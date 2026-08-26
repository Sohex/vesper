c eos.f
c KICO 04/08/08. Sub-routine for equation of state
c KICO 28/03/11. Updated to have two separate fuctions. eos calculates 
c density. eosd calculates vertical density gradient

c--------------------------------------------------
      subroutine eos(ec,t,s,z,ieos,rho)

      implicit none

      real ec(5), t, s, z, rho
      integer ieos 

      if(ieos.eq.0)then
c No thermobaricity term
        rho = ec(1)*t + ec(2)*s + ec(3)*t**2 + ec(4)*t**3
      elseif(ieos.eq.1)then
c Thermobaricity term is in
        rho = ec(1)*t + ec(2)*s + ec(3)*t**2 + ec(4)*t**3
     1            +ec(5)*t*z
      endif

      end

c---------------------------------------------------
c The two levels arrive as four SCALARS rather than as two two-element
c arrays. The only caller holds them in a four-dimensional array and would
c have to pass strided sections ts1(l,i,j,k:k+1) to fill an array dummy;
c gfortran packs each such section into a heap temporary before the call and
c frees it after, once per wet cell per level per ocean timestep. Scalars are
c passed by address and copy nothing. The arithmetic below is unchanged.
      subroutine eosd(ec,t1,t2,s1,s2,z,rdz,ieos,dzrho,tec)

      implicit none

      real ec(5), t1, t2, s1, s2, z, dzrho, tec
      integer ieos
      real tatw, rdz

c Calculate dzrho (vertical density gradient).

      tatw = 0.5*(t1 + t2)
      if(ieos.eq.0)then
c No thermobaricity term
        tec = - ec(1) - ec(3)*tatw*2 - ec(4)*tatw*tatw*3
      elseif(ieos.eq.1)then
c Thermobaricity term is in
        tec = - ec(1) - ec(3)*tatw*2 - ec(4)*tatw*tatw*3 - ec(5)*z
      endif
      dzrho = (ec(2)*(s2-s1) - tec*(t2-t1))*rdz

      end

