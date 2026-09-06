*
* subroutine velc for program goldstein variable depth
* fortran division replaced by reciprocals 190400
* first order boundary version but completely general for
* passages (for the first time) 10/7/1 and with simplified
* w calc assuming extended u array
* updated for generalised grid (RMA, 10/5/05)
*
c THREADED OVER j, in three worksharing loops rather than one.
c
c The three loops below are separate because the second and third read
c what the first wrote in a NEIGHBOURING row: the periodic copy takes
c u(1,imax,j,k) into u(1,0,j,k), and the vertical velocity reads
c u(2,i,j-1,k). The barrier at the end of each !$OMP DO is what orders
c them; merging them would be wrong at every row boundary between two
c threads.
c
c Inside each loop the algorithm is a column: the u and v recurrences run
c down k inside one (i,j) and nothing crosses to another column. rho,
c drag, dztau, dztav, rtv, rh and ub are read-only here, and u1 is read
c and written at the same index as u.
c
c dzu_col IS A LOCAL, and there is no COMMON member behind it. The two
c depth-integrated shear components are per-column scratch that this
c routine is the only reader and writer of anywhere in the tree. They
c were declared dzu(2,maxk) in ocean.cmn, which is one array for the
c whole process and so would have been shared between the threads; the
c declaration is gone from the header and from initialise_goldstein.
      subroutine velc

#include "ocean.cmn"

      real tv, tv1, tv2, tv4, tv5, sum(2)
      real dzu_col(2,maxk)
c      real tv3

      integer i, j, k, l

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP&   PRIVATE(i,j,k,l,tv,tv1,tv2,tv4,tv5,sum,dzu_col)

!$OMP DO SCHEDULE(STATIC)
      do 10 j=1,jmax
         do 20 i=1,imax
            do l=1,2
               sum(l) = 0
            enddo  
c u calc
            do 30 k=k1(i,j),kmax
               if(k1(i+1,j).gt.k)then
                  tv1 = 0
                  tv2 = 0
               else
                  tv2 = - (rho(i+1,j,k) - rho(i,j,k))*rdphi*rc(j)
                  if(max(k1(i,j-1),k1(i,j+1),k1(i+1,j-1),k1(i+1,j+1))
     1                     .le.k)then
                     tv1 = - c(j)*(rho(i+1,j+1,k) - rho(i+1,j-1,k)
     1                      + rho(i,j+1,k) - rho(i,j-1,k))*rds2(j)*0.25
                  elseif(max(k1(i,j-1),k1(i+1,j-1)).le.k)then
                     tv1 = - c(j)*(rho(i+1,j,k) - rho(i+1,j-1,k)
     1                      + rho(i,j,k) - rho(i,j-1,k))*rdsv(j-1)*0.5
                  elseif(max(k1(i,j+1),k1(i+1,j+1)).le.k)then
                     tv1 = - c(j)*(rho(i+1,j+1,k) - rho(i+1,j,k)
     1                      + rho(i,j+1,k) - rho(i,j,k))*rdsv(j)*0.5
                  else
                     tv1 = 0
                  endif
               endif
c
c v calc
c
               if(k1(i,j+1).gt.k)then
                  tv4 = 0
                  tv5 = 0
               else
                  tv4 = - cv(j)*(rho(i,j+1,k) - rho(i,j,k))*rdsv(j)
                  if(max(k1(i-1,j),k1(i-1,j+1),k1(i+1,j),k1(i+1,j+1))
     1                   .le.k)then
                     tv5 = - (rho(i+1,j+1,k) - rho(i-1,j+1,k)
     1                  + rho(i+1,j,k) - rho(i-1,j,k))*rdphi*0.25*rcv(j)
                  elseif(max(k1(i-1,j),k1(i-1,j+1)).le.k)then
                     tv5 = - (rho(i,j+1,k) - rho(i-1,j+1,k)
     1                  + rho(i,j,k) - rho(i-1,j,k))*rdphi*0.5*rcv(j)
                  elseif(max(k1(i+1,j),k1(i+1,j+1)).le.k)then
                     tv5 = - (rho(i+1,j+1,k) - rho(i,j+1,k)
     1                  + rho(i+1,j,k) - rho(i,j,k))*rdphi*0.5*rcv(j)
                  else
                     tv5 = 0
                  endif
               endif
c
c add wind
c
               if(k.eq.kmax)then
                  if(k1(i+1,j).le.k)then
                     tv1 = tv1 - dztau(2,i,j)
                     tv2 = tv2 - dztau(1,i,j)
                  endif
                  if(k1(i,j+1).le.k)then
                     tv4 = tv4 - dztav(2,i,j)
                     tv5 = tv5 - dztav(1,i,j)
                  endif
               endif
c
c integrate
c
               dzu_col(1,k) = - (s(j)*tv1 + drag(1,i,j)*tv2)*rtv(i,j)
               dzu_col(2,k) = - (drag(2,i,j)*tv4 - sv(j)*tv5)*rtv3(i,j)
               do l=1,2
                  if(k.eq.k1(i,j))then
                     u(l,i,j,k) = 0
                  else
                     u(l,i,j,k) = u(l,i,j,k-1) + dza(k-1)*(dzu_col(l,k)
     1                            + dzu_col(l,k-1))*0.5
                     sum(l) = sum(l) + dz(k)*u(l,i,j,k)
                  endif
               enddo  
   30       continue
            do k=k1(i,j),kmax
c
c add barotropic part and relax
c
               if(k1(i+1,j).le.k)then
                 u(1,i,j,k) = u(1,i,j,k) - sum(1)*rh(1,i,j) + ub(1,i,j)
                 u(1,i,j,k) = rel*u1(1,i,j,k) + (1.0 - rel)*u(1,i,j,k)
                 u1(1,i,j,k) = u(1,i,j,k)
               endif
               if(k1(i,j+1).le.k)then
                 u(2,i,j,k) = u(2,i,j,k) - sum(2)*rh(2,i,j) + ub(2,i,j)
                 u(2,i,j,k) = rel*u1(2,i,j,k) + (1.0 - rel)*u(2,i,j,k)
                 u1(2,i,j,k) = u(2,i,j,k)
               endif
            enddo
   20    continue
   10 continue
!$OMP END DO
c
c set boundary conditions 
c
!$OMP DO SCHEDULE(STATIC)
      do j=1,jmax
         do k=k1(1,j),kmax
            u(1,0,j,k) = u(1,imax,j,k)
         enddo
      enddo
!$OMP END DO
c
c calculate w
c
!$OMP DO SCHEDULE(STATIC)
      do j=1,jmax
         do i=1,imax
            tv = 0
            do k=k1(i,j),kmax-1
               tv1 = (u(1,i,j,k) - u(1,i-1,j,k))*rdphi*rc(j)
               tv2 = (u(2,i,j,k)*cv(j) - u(2,i,j-1,k)*cv(j-1))*rds(j)
               u(3,i,j,k) = tv           - dz(k)*(tv1 + tv2)
               tv = u(3,i,j,k)
c     if(k.lt.kmax.and.abs(u(3,i,j,k)).gt.1e-6)print*,i,j,k,u(3,i,j,k)
            enddo
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      end
