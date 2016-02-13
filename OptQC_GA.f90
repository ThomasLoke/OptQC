module csd_ga

use csd_tools
use m_mrgrnk
use rng

implicit none

type csd_ga_pool
	integer :: N, M
	integer :: pool_size, elite_size
	integer :: gen_lim									! Maximum number of attempts made to generate each chromosome - if fails just take the identity permutation
	double precision :: cross_rate, mut_rate 			! Real number between 0 and 1
	type(csd_solution_set), pointer :: css
	type(csd_generator), pointer :: cg
	integer :: obj_type
	double precision, allocatable :: X(:,:)
    double complex, allocatable :: Xc(:,:)
    integer :: ref_val									! Maximum number of gates - none of the costs should be greater than this value!!!
    double precision :: ft_total						! Total of fitness function values
	integer, allocatable :: chromosome_pool(:,:), temp_pool(:,:), chromosome_cost(:), chromosome_rank(:)
	double precision, allocatable :: chromosome_fitness(:)
	integer, allocatable :: Perm_temp(:)
contains
	procedure :: constructor => csd_ga_pool_constructor
	procedure :: destructor => csd_ga_pool_destructor
	procedure :: compute_cost => csd_ga_pool_compute_cost					! Can be used to recompute the optimal solution
	procedure :: fitness_function => csd_ga_pool_fitness_function
	procedure :: update_fitness => csd_ga_pool_update_fitness
	procedure :: evolve_rw => csd_ga_pool_evolve_rw
	procedure :: initialize_chromosomes => csd_ga_pool_initialize_chromosomes
	procedure :: sort_chromosomes => csd_ga_pool_sort_chromosomes
	procedure :: select_parents => csd_ga_pool_select_parents
	procedure :: perform_spc => csd_ga_pool_perform_spc
	procedure :: perform_crossover => csd_ga_pool_perform_crossover
	procedure :: perform_mutation => csd_ga_pool_perform_mutation
	procedure :: perform_step => csd_ga_pool_perform_step
	procedure :: print_state => csd_ga_pool_print_state
end type csd_ga_pool

contains

subroutine csd_ga_pool_constructor(this,css,cg,pool_size,elite_size,gen_lim,cross_rate,mut_rate)

implicit none
class(csd_ga_pool) :: this
type(csd_solution_set), target :: css
type(csd_generator), target :: cg
integer :: pool_size, elite_size, gen_lim
double precision :: cross_rate, mut_rate

integer :: M

M = css%M
this%N = css%N
this%M = M
this%pool_size = pool_size
this%elite_size = elite_size
this%gen_lim = gen_lim
this%cross_rate = cross_rate
this%mut_rate = mut_rate
this%css => css
this%cg => cg
this%obj_type = css%arr(3)%obj_type
if(this%obj_type == 0) then
	allocate(this%X(M,M))
	this%X = css%arr(3)%X
else
	allocate(this%Xc(M,M))
	this%Xc = css%arr(3)%Xc
end if
this%ref_val = css%csdr_ss_ct
allocate(this%chromosome_pool(pool_size,M))
allocate(this%temp_pool(pool_size,M))
allocate(this%chromosome_cost(pool_size))
allocate(this%chromosome_rank(pool_size))
allocate(this%chromosome_fitness(pool_size))
allocate(this%Perm_temp(M))

end subroutine csd_ga_pool_constructor

subroutine csd_ga_pool_destructor(this)

implicit none
class(csd_ga_pool) :: this

this%N = 0
this%M = 0
this%pool_size = 0
this%elite_size = 0
this%gen_lim = 0
this%cross_rate = 0.0d0
this%mut_rate = 0.0d0
nullify(this%css)
nullify(this%cg)
if(this%obj_type == 0) then
	deallocate(this%X)
else
	deallocate(this%Xc)
end if
this%ref_val = 0
deallocate(this%chromosome_pool)
deallocate(this%temp_pool)
deallocate(this%chromosome_cost)
deallocate(this%chromosome_rank)
deallocate(this%chromosome_fitness)
deallocate(this%Perm_temp)

end subroutine csd_ga_pool_destructor

function csd_ga_pool_compute_cost(this,Perm)

implicit none
class(csd_ga_pool) :: this
integer :: Perm(:)
integer :: csd_ga_pool_compute_cost

! Construct the P matrices
call permlisttomatrixtr(this%M,Perm,this%css%arr(2)%X)
call permlisttomatrix(this%M,Perm,this%css%arr(4)%X)
! Update the U' matrix - done with reference to the stored copy
if(this%css%arr(3)%obj_type == 0) then
    call ApplyPerm(this%M,this%X,this%css%arr(3)%X,Perm)
else
    call ApplyPerm_CPLX(this%M,this%Xc,this%css%arr(3)%Xc,Perm)
end if
! Calculate the new cost
call this%css%run_csdr(this%cg)
csd_ga_pool_compute_cost = this%css%csdr_ss_ct

end function csd_ga_pool_compute_cost

subroutine csd_ga_pool_fitness_function(this)

implicit none
class(csd_ga_pool) :: this

this%chromosome_fitness = maxval(this%chromosome_cost) - this%chromosome_cost + 1.0d0

end subroutine csd_ga_pool_fitness_function

subroutine csd_ga_pool_update_fitness(this)

implicit none
class(csd_ga_pool) :: this

integer :: i

! Only update cost values for non-elite chromosomes
do i = this%elite_size+1, this%pool_size
	this%chromosome_cost(i) = this%compute_cost(this%chromosome_pool(i,:))
end do
call this%fitness_function()

end subroutine csd_ga_pool_update_fitness

subroutine csd_ga_pool_evolve_rw(this,Perm,cost)

implicit none
class(csd_ga_pool) :: this
integer :: Perm(:), cost

integer :: M, i, cost_temp, s1, s2, temp

M = this%M
do i = 1, this%gen_lim
	! Duplicate the current solution
	this%Perm_temp = Perm
	cost_temp = cost
	! Mutate the old solution to produce a new solution
	s1 = rng_inst%rint(M)
	s2 = rng_inst%rint(M)
	do while (s2 == s1)
		s2 = rng_inst%rint(M)
	end do
	temp = this%Perm_temp(s1)
	this%Perm_temp(s1) = this%Perm_temp(s2)
	this%Perm_temp(s2) = temp
	! Compute the cost of the new solution
	cost_temp = this%compute_cost(this%Perm_temp)
	! If the new solution is more than or at least as optimal as the old solution, accept the new solution
	if(cost_temp <= cost) then
		Perm = this%Perm_temp
		cost = cost_temp
	end if
end do

end subroutine csd_ga_pool_evolve_rw

subroutine csd_ga_pool_initialize_chromosomes(this)

implicit none
class(csd_ga_pool) :: this

integer :: M, ref_val, i

M = this%M
ref_val = this%ref_val
! Copy the identity permutation to all the chromosomes
do i = 1, M
	this%chromosome_pool(:,i) = i
end do
this%chromosome_cost(:) = ref_val
! Fix the first chromosome as the identity permutation, but generate the rest as normal
do i = 2, this%pool_size
	call this%evolve_rw(this%chromosome_pool(i,:),this%chromosome_cost(i))
end do
! Finally, compute the corresponding fitness values
call this%fitness_function()
! Sort the corresponding results
call this%sort_chromosomes()

end subroutine csd_ga_pool_initialize_chromosomes

! subroutine csd_ga_pool_initialize_chromosomes(this)

! implicit none
! class(csd_ga_pool) :: this

! integer :: M, ref_val, i, j, temp
! logical :: found

! M = this%M
! ref_val = this%ref_val
! ! Fix the first chromosome as the identity permutation
! do i = 1, M
! 	this%chromosome_pool(1,i) = i
! end do
! this%chromosome_cost(1) = ref_val
! ! Generate the rest as normal
! do i = 1, this%pool_size
! 	found = .false.
! 	do j = 1, this%gen_lim
! 		call perm_generate(M,this%chromosome_pool(i,:))
! 		temp = this%compute_cost(this%chromosome_pool(i,:))
! 		if(.true. .or. temp <= ref_val) then
! 			this%chromosome_cost(i) = temp
! 			found = .true.
! 			exit
! 		end if
! 	end do
! 	! If cannot generate a better or equal chromosome, take the identity permutation
! 	if(found == .false.) then
! 		this%chromosome_pool(i,:) = this%chromosome_pool(1,:)
! 		this%chromosome_cost(i) = ref_val
! 	end if
! end do
! ! Finally, compute the corresponding fitness values
! call this%fitness_function()
! ! Sort the corresponding results
! call this%sort_chromosomes()

! end subroutine csd_ga_pool_initialize_chromosomes

subroutine csd_ga_pool_sort_chromosomes(this)

implicit none
class(csd_ga_pool) :: this

! Compute the ranking of each chromosome - fittest (i.e. smallest cost and highest fitness) is rank 1
call mrgrnk(this%chromosome_cost,this%chromosome_rank)
! Sort the population from fittest to least fit
this%chromosome_pool = this%chromosome_pool(this%chromosome_rank,:)
! Sort the associated cost and fitness as well
this%chromosome_cost = this%chromosome_cost(this%chromosome_rank)
this%chromosome_fitness = this%chromosome_fitness(this%chromosome_rank)

end subroutine csd_ga_pool_sort_chromosomes

subroutine csd_ga_pool_select_parents(this,parents)

implicit none
class(csd_ga_pool) :: this
integer :: parents(2)

double precision :: ft_r(2), ft, spacing
integer :: i, pos

! Assume ft_total has been calculated beforehand
! Use stochastic universal sampling with 2 points - note that inbreeding can occur
spacing = 0.5d0 * this%ft_total
ft_r(1) = rng_inst%r01() * spacing
ft_r(2) = ft_r(1) + spacing
do i = 1, 2
	ft = 0.0d0
	pos = 0
	do while(ft < ft_r(i))
		pos = pos+1
		ft = ft + this%chromosome_fitness(pos)
	end do
	if(pos == 0) pos = 1		! Just in case ft_r = 0
	parents(i) = pos
end do

end subroutine csd_ga_pool_select_parents

subroutine csd_ga_pool_perform_spc(this,c_parents,c_new)

implicit none
class(csd_ga_pool) :: this
integer :: c_parents(:,:), c_new(:)

integer :: i, j, idx_new, r, ref
logical :: found

! Select the crossover point
r = rng_inst%rint(this%M)
! Copy over first parent up to the crossover point
c_new(1:r) = c_parents(1,1:r)
! The rest is taken in the same order as in the second parent
idx_new = r+1
do i = 1, this%M
	ref = c_parents(2,i)
	! After this loop, if found = .true., then a match was found
	found = .false.
	do j = 1, r
		if(c_new(j) == ref) then
			found = .true.
			exit
		end if
	end do
	! If no match found, add it to the new chromosome
	if(found == .false.) then
		c_new(idx_new) = ref
		idx_new = idx_new + 1
		! If new chromosome has been completely filled, procedure is complete, so skip the rest
		if(idx_new > this%M) then
			exit
		end if
	end if
end do

end subroutine csd_ga_pool_perform_spc

subroutine csd_ga_pool_perform_crossover(this)

implicit none
class(csd_ga_pool) :: this

integer :: i, parents(2)

! Assume sorted beforehand
! Clear the temporary pool - use to construct the new chromosome pool
this%temp_pool = 0
! Copy the elite chromosomes over first - retain cost values (post-sorting)
this%temp_pool(1:this%elite_size,:) = this%chromosome_pool(1:this%elite_size,:)
! Compute total fitness of the population
this%ft_total = sum(this%chromosome_fitness)
! Perform crossovers between selected chromosomes to form the rest of the chromosome pool
do i = this%elite_size+1, this%pool_size
	! Select parent chromosomes
	call this%select_parents(parents)
	! Perform single point crossover with the parent chromosomes to form the new chromosome
	call this%perform_spc(this%chromosome_pool(parents,:),this%temp_pool(i,:))
end do
! After new chromosome pool has been formed, copy it back from the temp pool
this%chromosome_pool = this%temp_pool

end subroutine csd_ga_pool_perform_crossover

subroutine csd_ga_pool_perform_mutation(this)

implicit none
class(csd_ga_pool) :: this

integer :: M, i, s1, s2, temp

M = this%M
! Note: No mutation is performed on the elite chromosomes
do i = this%elite_size+1, this%pool_size
	if(rng_inst%r01() < this%mut_rate) then
		! Perform normal random mutation - interchange two random positions in the chromosome
		s1 = rng_inst%rint(M)
		s2 = rng_inst%rint(M)
		do while (s2 == s1)
		    s2 = rng_inst%rint(M)
		end do
		temp = this%chromosome_pool(i,s1)
		this%chromosome_pool(i,s1) = this%chromosome_pool(i,s2)
		this%chromosome_pool(i,s2) = temp
		! No need to update cost values
	end if
end do	

end subroutine csd_ga_pool_perform_mutation

subroutine csd_ga_pool_perform_step(this)

implicit none
class(csd_ga_pool) :: this

call this%perform_crossover()
call this%perform_mutation()
call this%update_fitness()
call this%sort_chromosomes()

end subroutine csd_ga_pool_perform_step

subroutine csd_ga_pool_print_state(this)

implicit none
class(csd_ga_pool) :: this

integer :: i

do i = 1, this%pool_size
	write(*,'(a,i4,a,i8,a,f15.9)')'Chromosome ',i,': ',this%chromosome_cost(i),' / ',this%chromosome_fitness(i)
end do

end subroutine csd_ga_pool_print_state

end module csd_ga