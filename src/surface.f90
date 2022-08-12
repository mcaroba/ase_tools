module surface_tools

  contains

  subroutine surface_atoms(positions, cell, min_r, max_r, n_tries, cluster, is_surface)

    implicit none

    real*8, intent(in) :: positions(:,:), min_r, max_r, cell(1:3,1:3)
    integer, intent(in) :: n_tries
    logical, intent(in) :: cluster
    logical, intent(out), dimension(1:size(positions,1)) :: is_surface

    real*8 :: rmin(1:3), rmax(1:3), d, dv(1:3), a_box(1:3), b_box(1:3), c_box(1:3)
    real*8, allocatable :: rand(:,:)
    logical, allocatable :: too_close(:)
    integer :: n_atoms, i, j

    n_atoms = size(positions, 1)

    a_box(1:3) = cell(1,1:3)
    b_box(1:3) = cell(2,1:3)
    c_box(1:3) = cell(3,1:3)

    allocate( rand(1:n_tries, 1:3) )
    call random_number(rand)

    allocate( too_close(1:n_tries) )
    too_close = .false.

    if( cluster )then
!     If the structure is a cluster we disregard the lattice vectors and use a reduced orthogonal
!     box for the MC sampling
      rmin = positions(1, 1:3)
      do i = 2, n_atoms
        do j = 1, 3
          if( positions(i, j) < rmin(j) )then
            rmin(j) = positions(i, j)
          end if
        end do
      end do
      rmax = positions(1, 1:3)
      do i = 2, n_atoms
        do j = 1, 3
          if( positions(i, j) > rmax(j) )then
            rmax(j) = positions(i, j)
          end if
        end do
      end do

      rmin = rmin - min_r
      rmax = rmax + min_r

      do i = 1, n_tries
        rand(i,1:3) = rmin(1:3) + rand(i,1:3)*(rmax(1:3)-rmin(1:3))
      end do

      do i = 1, n_tries
        do j = 1, n_atoms
          dv = positions(j,:) - rand(i,:)
          d = sqrt(dot_product(dv,dv))
          if( d < min_r )then
            too_close(i) = .true.
            exit
          end if
        end do
      end do

      is_surface = .false.
      do j = 1, n_atoms
        do i = 1, n_tries
          if( .not. too_close(i) )then
            dv = positions(j,:) - rand(i,:)
            d = sqrt(dot_product(dv,dv))
            if( d < max_r )then
              is_surface(j) = .true.
              exit
            end if
          end if
        end do
      end do
    else
      write(*,*) "Not implemented!"
      stop
      do i = 1, n_tries
        rand(i,1:3) = rand(i,1)*a_box(1:3) + rand(i,2)*b_box(1:3) + rand(i,3)*c_box(1:3)
      end do
    end if


    deallocate( rand, too_close )

  end subroutine

end module
