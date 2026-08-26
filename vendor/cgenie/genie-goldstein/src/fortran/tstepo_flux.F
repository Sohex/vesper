c **********************************************************************
c *** ocean flux calculation *******************************************
c **********************************************************************
c
c THE PARALLEL AXIS IS j, AND WHY IT IS j RATHER THAN k OR i.
c
c This routine is half the ocean's retired instructions, so it is the
c routine the threading is for. Its nest is k outside, then j, then i, and
c each index carries a different kind of dependence:
c
c   k  fb(l,i,j) holds the flux through the BOTTOM face and is filled at
c      level k for use at level k+1. That is a real dependence over the
c      whole domain, so k stays sequential and the levels are separated by
c      the barrier at the end of the !$OMP DO.
c
c   i  fw(l) holds the flux through the WESTERN face and is filled by the
c      previous cell in the row. That is contained within one j, so it does
c      not constrain a decomposition over j at all.
c
c   j  fs(l,i) holds the flux through the SOUTHERN face and was filled by
c      the previous ROW. That is the one dependence a decomposition over j
c      has to remove, and it is the recoverable kind: the northern flux of
c      row j-1 is a pure function of u, ts1 and k1 at that row, all of
c      which are read-only for the whole routine. So each row RECOMPUTES
c      its own southern face rather than inheriting it. That costs one
c      extra meridional flux per cell and removes the dependence exactly:
c      the expression below is the same expression, at j-1, that the
c      recycling filled fs from.
c
c ts1 is not written anywhere in this routine -- tstepo updates it after
c this returns -- which is what makes that recomputation exact rather than
c approximate.
c
c THE TWO REDUCTIONS ARE BOTH ORDER-INDEPENDENT. limps counts flux-limited
c points and dmax takes a maximum, so neither depends on the thread count
c and the result is reproducible rather than merely correct. They are
c accumulated into locals because limps and dmax themselves are members of
c a COMMON block.
c
c diffv IS INITIALISED INSIDE THE PARALLEL REGION and that is not
c cosmetic. It is assigned per cell only when iediff is 1 or 2 AND the
c cell is wet and below the surface; everywhere else it keeps the value it
c had, which for iediff = 0 is diff(2) from before the loop. A private
c copy that inherited nothing would be whatever was on the stack.

      subroutine tstepo_flux

#include "ocean.cmn"

      real tv, ups(3), ups0, pec(3)

      real fe(maxl), fw(maxl), fn(maxl), fs(maxl,maxi), fa(maxl)
     1   , fb(maxl,maxi,maxj), fwsave(maxl)
      integer i, j, k, l
      parameter(ups0 = 0.)

c KICO ediff calc needs
      real diffv

      real tec, scc, dzrho, rdzrho, slim, tv1
      real dxrho(4), dxts(maxl,4), dyrho(4), dyts(maxl,4), dzts(maxl)
      integer ina, nnp, knp
      integer loc_lmax
c reduction accumulators, see the header
      integer loc_limps
      real loc_dmax

c ------------------------------------------------------------ c
c INITIALIZE VARIABLES
c ------------------------------------------------------------ c
      scc = 0.0
c ------------------------------------------------------------ c

      if(diso)then
         scc = ec(2)
         limps = 0
      endif

      dmax = 0
      loc_limps = 0
      loc_dmax = 0

* 2nd order explicit step

c lower boundary fluxes
      do 220 j=1,jmax
         do 221 i=1,imax
            do 222 l=1,lmax
               fb(l,i,j) = 0
 222        continue
 221     continue
 220  continue

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&   PRIVATE(i,j,k,l,tv,ups,pec,fe,fw,fn,fs,fa,fwsave,diffv,
!$OMP&           tec,dzrho,rdzrho,slim,tv1,dxrho,dxts,dyrho,dyts,dzts,
!$OMP&           ina,nnp,knp,loc_lmax)
!$OMP&   REDUCTION(+:loc_limps) REDUCTION(MAX:loc_dmax)

c KICO initialise diffv
      diffv=diff(2)
      rdzrho = 0.0

      do 100 k=1,kmax
!$OMP DO SCHEDULE(STATIC)
         do 101 j=1,jmax
c southern boundary fluxes: the northern flux of row j-1, RECOMPUTED
c rather than inherited from the previous iteration, which is what lets
c this loop be divided over threads. See the header.
            if(j.eq.1)then
               do i=1,imax
                  do l=1,lmax
                     fs(l,i) = 0
                  enddo
               enddo
            else
               do i=1,imax
                  pec(2) = u(2,i,j-1,k)*dsv(min(j-1,jmax-1))/diff(1)
                  ups(2) = pec(2) / (2.0 + abs(pec(2)))
                  do l=1,lmax
                     if(k.lt.max(k1(i,j-1),k1(i,j)))then
                        fs(l,i) = 0
                     else
                        fs(l,i) = cv(j-1)*u(2,i,j-1,k)*((1.-ups(2))
     &                   *ts1(l,i,j,k) + (1.+ups(2))*ts1(l,i,j-1,k))*0.5
                        fs(l,i) = fs(l,i) - cv2(j-1)*(ts1(l,i,j,k) -
     1                       ts1(l,i,j-1,k))*diff(1)
                     endif
                  enddo
               enddo
            endif
c western boundary fluxes
            i = 1
            pec(1) = u(1,imax,j,k)*dphi/diff(1)
            ups(1) = pec(1) / (2.0 + abs(pec(1)))
            do 210 l=1,lmax
               if (k.ge.max(k1(imax,j),k1(1,j)))then
c western doorway
                  fw(l) = u(1,imax,j,k)*rc(j)*((1.-ups(1))*ts1(l,1,j,k)
     1                  + (1.+ups(1))*ts1(l,imax,j,k))*0.5
                  fw(l) = fw(l) - (ts1(l,1,j,k) - ts1(l,imax,j,k))
     1                       *rc2(j)*diff(1)
               else
                  fw(l) = 0
               endif
               fwsave(l) = fw(l)
  210       continue
            do 102 i=1,imax

c KICO calculate local vertical diffusivity
              if(k.ge.k1(i,j).and.k.lt.kmax)then
c First get vertical density gradient (also needed for isopycnal diff)
                call eosd(ec,ts1(1,i,j,k),ts1(1,i,j,k+1),
     1             ts1(2,i,j,k),ts1(2,i,j,k+1),zw(k),
     2             rdza(k),ieos,dzrho,tec)
                if(dzrho.lt.-1e-12)then
                  rdzrho = 1.0/dzrho
                else
                  rdzrho = -1e12
                endif
                dzrho_test(i,j,k) = dzrho
                if((iediff.gt.0).AND.(iediff.lt.3))then
c Value of diffv fine for applying diffusivity, but peclet number calc
c is a 1st order approximation if diffusivity is variable.
                   if(ediffpow2i.eq.0)then
                     diffv=ediff0+ediff1(i,j,k)
                   elseif(ediffpow2i.eq.1)then
                     diffv=ediff0+ediff1(i,j,k)*(-rdzrho)
                   elseif(ediffpow2i.eq.2)then
                     diffv=ediff0+ediff1(i,j,k)*sqrt(-rdzrho)
                   else
                     diffv=ediff0+ediff1(i,j,k)*((-rdzrho)**ediffpow2)
                   endif
                   if(diffv.gt.diffmax(k+1)) diffv=diffmax(k+1)
                   diffv_test(i,j,k) = diffv
                 endif
               endif
               pec(1) = u(1,i,j,k)*dphi/diff(1)
               ups(1) = pec(1) / (2.0 + abs(pec(1)))
c rather untidy mask to avoid undefined dsv at jmax nre
               pec(2) = u(2,i,j,k)*dsv(min(j,jmax-1))/diff(1)
               ups(2) = pec(2) / (2.0 + abs(pec(2)))
               pec(3) = u(3,i,j,k)*dza(k)/diffv
               ups(3) = pec(3) / (2.0 + abs(pec(3)))
               do l=1,lmax
c flux to east
                  if(i.eq.imax)then
c eastern edge(doorway or wall)
                     fe(l) = fwsave(l)
                  elseif(k.lt.max(k1(i,j),k1(i+1,j)))then
                     fe(l) = 0
                  else
                     fe(l) = u(1,i,j,k)*rc(j)*((1.-ups(1))
     &                   *ts1(l,i+1,j,k) + (1.+ups(1))*ts1(l,i,j,k))*0.5
                     fe(l) = fe(l) - (ts1(l,i+1,j,k) - ts1(l,i,j,k))
     1                       *rc2(j)*diff(1)
                  endif
c flux to north
                  if(k.lt.max(k1(i,j),k1(i,j+1)))then
                     fn(l) = 0
                  else
                     fn(l) = cv(j)*u(2,i,j,k)*((1.-ups(2))
     &                   *ts1(l,i,j+1,k) + (1.+ups(2))*ts1(l,i,j,k))*0.5
                     fn(l) = fn(l) - cv2(j)*(ts1(l,i,j+1,k) -
     1                       ts1(l,i,j,k))*diff(1)
                  endif
c flux above
                  if(k.lt.k1(i,j))then
                     fa(l) = 0
                  elseif(k.eq.kmax)then
                     fa(l) = ts(l,i,j,kmax+1) 
                  else
                     fa(l) = u(3,i,j,k)*((1.-ups(3))*ts1(l,i,j,k+1) +
     1                       (1.+ups(3))*ts1(l,i,j,k))*0.5
                     fa(l) = fa(l) - (ts1(l,i,j,k+1) - ts1(l,i,j,k))
     1                       *rdza(k)*diffv
                  endif
               enddo
c     isoneutral diffusion
               if(diso)then
c limit isoneutral diffusion to only T and S?
                  if (ctrl_diso_TSonly) then
                     loc_lmax = 2
                  else
                     loc_lmax = lmax                     
                  end if
c check for valid ocean level (k)
                  if(k.ge.k1(i,j).and.k.lt.kmax)then
                     if(dzrho.lt.-1e-12)then
                        tv1 = 0.0
c tracer loop          
                        do knp=0,1
                           do nnp=0,1
                              ina = 1+nnp + 2*knp
c phi derivatives
                              do l=1,loc_lmax
                                 if(k+knp.ge.k1(i-1+2*nnp,j))then
                                    dxts(l,ina) = (ts1(l,i+nnp,j,k+knp)
     2                                   - ts1(l,i+nnp-1,j,k+knp))
     3                                   *rc(j)*rdphi
                                 else
                                    dxts(l,ina) = 0.
                                 endif
c s-derivatives
                                 if(k+knp.ge.k1(i,j-1+2*nnp))then
                                    dyts(l,ina) = (ts1(l,i,j+nnp,k+knp)
     2                                   - ts1(l,i,j+nnp-1,k+knp))
     &                                   *cv(j-1+nnp)*rdsv(j+nnp-1)
                                 else
                                    dyts(l,ina) = 0.
                                 endif
                              enddo
                              dxrho(ina) = scc*dxts(2,ina)
     1                             -tec*dxts(1,ina)
                              dyrho(ina) = scc*dyts(2,ina)
     1                             -tec*dyts(1,ina)
c calculate diagonal part
                              tv1 = tv1 + dxrho(ina)*dxrho(ina)
     1                             + dyrho(ina)*dyrho(ina)
                           enddo
                        enddo
                        tv1 = 0.25*tv1*rdzrho*rdzrho
c limit flux by factor slim for large slope
                        if(tv1.gt.ssmax(k))then
                           slim = ssmax(k)*ssmax(k)/(tv1*tv1)
c count flux-limited points
                           loc_limps = loc_limps + 1
                        else
                           slim = 1.0
                        endif
                        tv1 = tv1*slim*diff(1)*rdza(k)
c test vertical diffusion number
                        tv = tv1*dt(k)*rdza(k)
                        if(tv.gt.loc_dmax)then
                           loc_dmax = tv
                        endif
                        do l=1,loc_lmax
                           dzts(l) = (ts1(l,i,j,k+1) 
     1                          - ts1(l,i,j,k))*rdza(k)
c add isoneutral vertical flux
                           tv = 0
                           do ina=1,4
                              tv = tv + (2*dzrho*dxts(l,ina)
     1                             - dxrho(ina)*dzts(l))*dxrho(ina) 
     2                             + (2*dzrho*dyts(l,ina)
     3                             - dyrho(ina)*dzts(l))*dyrho(ina)
                           enddo
                           tv = 0.25*slim*diff(1)*tv/(dzrho*dzrho)
                           fa(l) = fa(l) + tv
                        enddo
                     endif
                  endif
               endif
               do l=1,lmax
                  tv = 0
                  if(k.ge.k1(i,j))then
                     ts(l,i,j,k) = ts1(l,i,j,k) - dt(k)*( - tv +
     1                           (fe(l) - fw(l))*rdphi
     2                           + (fn(l) - fs(l,i))*rds(j)
     3                           + (fa(l) - fb(l,i,j))*rdz(k))
                  endif
                  fw(l) = fe(l)
c fs(l,i) = fn(l) was here. It is a dead store now: row j+1 recomputes
c its own southern face from ts1 and u rather than inheriting this one,
c which is what removes the j dependence.
                  fb(l,i,j) = fa(l)
               enddo

               call eos(ec,ts(1,i,j,k),ts(2,i,j,k),zro(k),
     1            ieos,rho(i,j,k))
 102        continue
 101     continue
!$OMP END DO
 100  continue
!$OMP END PARALLEL

      limps = limps + loc_limps
      dmax = loc_dmax

      end
