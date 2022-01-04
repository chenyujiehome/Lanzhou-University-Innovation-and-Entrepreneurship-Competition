module ran_mod   
! module contains four functions   
! ran1 and ran2 return a uniform random number between 0-1   
! normal1 and normal 2 return a normal distribution   
contains  
   
    function ran1()    !returns random number between 0 - 1   
    implicit none  
    integer :: flag 
    double precision :: ran1 
    save flag 
    data flag /0/  
 
      if(flag==0) then 
        call random_seed(); flag = 1 
      endif 
      call random_number(ran1)                   ! built in fortran 90 random number function   
 
    return 
    end function ran1  
 
    function ran2(idum) 
    implicit none 
    integer, parameter :: IM1=2147483563, IM2=2147483399 
    integer, parameter :: IMM1=IM1-1, IA1=40014, IA2=40692 
    integer, parameter :: IQ1=53668, IQ2=52774, IR1=12211, IR2=3791 
    integer, parameter :: NTAB=32, NDIV=1+IMM1/NTAB 
    integer, dimension (NTAB) :: iv 
    integer :: idum, idum2, iy, i, j, k 
    real, parameter :: EPS=1.2D-7, RNMX=1.d0-EPS, AM=1.d0/IM1 
    double precision :: ran2, RNMK 
    save iv, iy, idum2 
    data idum2/123456789/, iv/NTAB*0/, iy/0/ 
 
      if(idum.le.0)then 
        idum = max(-idum,1) 
        idum2 = idum 
        do j = NTAB+8,1,-1 
          k = idum/IQ1 
          idum = IA1*(idum-k*IQ1)-k*IR1 
          if(idum.lt.0) idum = idum + IM1 
          if(j.le.NTAB) iv(j) = idum 
        end do 
        iy = iv(1) 
      end if 
      k = idum/IQ1 
      idum = IA1*(idum-k*IQ1)-k*IR1 
      if(idum.lt.0) idum = idum + IM1 
      k = idum2/IQ2 
      idum2 = IA2*(idum2-k*IQ2)-k*IR2 
            if(idum2.lt.0) idum2 = idum2 + IM2 
      j = 1+iy/NDIV 
      iy = iv(j)-idum2 
      iv(j) = idum 
            if(iy.lt.1) iy = iy + IMM1 
      ran2 = min(AM*iy,RNMX) 
 
    return 
    end function 
 
    function normal1(mean,sigma) !returns a normal distribution   
    implicit none  
    integer :: flag 
    double precision :: normal1, tmp, mean, sigma   
    double precision :: fac, gsave, rsq, r1, r2   
    save flag, gsave  
    data flag /0/  
 
      if (flag.eq.0) then  
        rsq=2.0d0  
        do while(rsq.ge.1.0d0.or.rsq.eq.0.0d0) ! new from for do  
          r1 = 2.0d0 * ran1() - 1.0d0  
          r2 = 2.0d0 * ran1() - 1.0d0  
                    rsq = r1 * r1 + r2 * r2   
        enddo  
        fac = sqrt(-2.0d0*log(rsq)/rsq)  
        gsave = r1 * fac  
        tmp = r2 * fac  
        flag = 1   
      else  
        tmp = gsave  
        flag = 0  
      endif  
      normal1 = tmp * sigma + mean   
 
    return  
    end function normal1  
 
    function normal2(mean,sigma) 
    implicit none 
    integer :: flag 
    double precision, parameter :: pi = 3.141592653589793239   
    double precision :: u1, u2, y1, y2, normal2, mean, sigma 
    save flag 
    data flag /0/ 
 
      u1 = ran1(); u2 = ran1() 
      if (flag.eq.0) then 
        y1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2) 
        normal2 = mean + sigma*y1 
        flag = 1 
      else 
        y2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2) 
        normal2 = mean + sigma*y2 
        flag = 0 
      endif 
 
    return 
    end function normal2 
 
end module ran_mod 
 
program main
use ran_mod
implicit none  
real*8 :: mean,sigma,ur1,ur2,nr1,nr2
integer :: i,iseed

  iseed = 12345; mean = 0.d0; sigma = 1.d0
  open(66,file='d:\test\svd3\Rdata.txt')
  do i = 1, 100000
 !   ur1 = ran1(); ur2 = ran2(iseed)
 !   nr1 = normal1(mean,sigma)
    nr2 = normal2(mean,sigma)
 !   write(66,'(4(1X,F9.6))') ur1,ur2,nr1,nr2
    write(66,'(1X,F9.6)') nr2
  enddo
  close(66)

endprogram main
