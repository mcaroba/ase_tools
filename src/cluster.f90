module cluster_module

  implicit none

  contains

    recursive subroutine find_neighbors(i, j, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)

      implicit none

!     Input variables
      integer, intent(in) :: i, j, n_atoms, cluster
      logical, intent(in) :: bonded(:,:)
!     Inout variables
      integer, intent(inout) :: atom_belongs_to_cluster(:)
      logical, intent(inout) :: atom_visited(:)
!     Internal variables
      integer :: k

      if( .not. atom_visited(j) .and. i /= j .and. bonded(i,j) )then
        atom_visited(j) = .true.
        atom_belongs_to_cluster(j) = cluster
        do k = 1, n_atoms
          call find_neighbors(j, k, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)
        end do
      end if
    end subroutine


    subroutine get_distance(pos1, pos2, lx, ly, lz, d)

!   Returns distance d between atoms i and j according
!   to minimum image convention

      implicit none

      real*8, intent(out) :: d
      real*8, intent(in) :: lx, ly, lz, pos1(1:3), pos2(1:3)
      real*8 :: x, y, z

      x = pos2(1) - pos1(1)
      y = pos2(2) - pos1(2)
      z = pos2(3) - pos1(3)

      if( dabs(x) > lx/2.d0)then
        x = x - sign(lx, x)
      end if
      if( dabs(y) > ly/2.d0)then
        y = y - sign(ly, y)
      end if
      if( dabs(z) > lz/2.d0)then
        z = z - sign(lz, z)
      end if

      d = sqrt( x**2 + y**2 + z**2)

    end subroutine



    subroutine cluster_atoms(positions, lx, ly, lz, bonding_cutoff, atom_belongs_to_cluster)

      implicit none

!     Input variables
      real*8, intent(in) :: positions(:,:), lx, ly, lz, bonding_cutoff(:)
!     Output variables
      integer, intent(out), dimension(1:size(positions,2)) :: atom_belongs_to_cluster
!     Internal variables
      real*8 :: d
      integer :: n_atoms, cluster, i, j
      logical, allocatable :: bonded(:,:), atom_visited(:)

      n_atoms = size(positions,2)

      allocate( bonded(1:n_atoms, 1:n_atoms) )
      bonded = .false.

      !$omp parallel do private(i,j,d)
      do i = 1, n_atoms
        do j = i+1, n_atoms
          call get_distance(positions(1:3, i), positions(1:3, j), lx, ly, lz, d)
          if( d < (bonding_cutoff(i)+bonding_cutoff(j))/2.d0 )then
            bonded(i, j) = .true.
            bonded(j, i) = .true.
          end if
        end do
      end do

      atom_belongs_to_cluster = 0
      allocate( atom_visited(1:n_atoms) )
      atom_visited = .false.
      cluster = 0
      do i = 1, n_atoms
!     Find all the atoms that can be connected to i and put them on the same cluster
        if( .not. atom_visited(i) )then
          cluster = cluster + 1
          atom_belongs_to_cluster(i) = cluster
          atom_visited(i) = .true.
          do j = 1, n_atoms
            call find_neighbors(i, j, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)
          end do
        end if
      end do

      deallocate( bonded, atom_visited )

    end subroutine


end module
