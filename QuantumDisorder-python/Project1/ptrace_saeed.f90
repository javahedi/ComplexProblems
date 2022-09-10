subroutine get_trace(Psi, da, db, rho)  ! Returns the right partial trace (over b), for a bi-partite matrix
implicit none
integer, intent(in) :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8), intent(in) :: Psi(1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8), intent(out) :: rho(1:da,1:da)  !  Reduced matrix
integer :: i,j,m  ! Auxiliary variables for counters

rho = 0.d0
do i = 1, da    
    do j = i, da
          do m = 1, db    
              rho(i,j) = rho(i,j) + Psi(m+i*db)*Psi(m+i*db+(j-i)*db) 
          enddo
          if ( j /= i ) rho(j,i) = conjg(rho(i,j))
    enddo 
enddo

end

