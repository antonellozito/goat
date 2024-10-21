!======================================================================!
!                                                                      !
!                               DOCUMENTATION                          !
!                                                                      !
!======================================================================!

! Description
!============
! This module contains procedures to solve quadratic problems (QP). 
! Linear (in)equality constraints are allowed. In the case of inequality
! constraints, the solution will proceed iteratively in order to 
! determine the active set. This module uses the hessian approximation
! module to represent the hessian. Gradients are assumed to be given as
! a one-dimensional array. 

! Important
!==========
! The current implementation is naive and non-optimized. 
! It should only be used for problems that are relatively
! small in the sense that they don't require (massive) parallellization
! or (very) special solver techniques. Some large scale problems may 
! be handled quite effectively, depending on the hessian representation
! and solution procedure (e.g. unconstrained problems with 
! limited-memory bfgs). 

module optmod_qp 

    ! Initialize
    !============
    ! Load modules
    use mod_precision 
    use mod_sparseinterface
    use mod_linearsolverinterface
    use optmod_hessianapproximation

    ! The usual
    implicit none
    save
    private 

    ! Set general solver routines to public
    public SolveQPDirect

    !==================================================================!
    !                                                                  !
    !                               TYPES                              !
    !                                                                  !
    !==================================================================!

    !==================================================================!
    !                                                                  !
    !                            INTERFACES                            !
    !                                                                  !
    !==================================================================!

    ! Direct QP solver interface
    interface SolveQPDirect 
        module procedure SolveQPDirectNocon
        module procedure SolveQPDirectEqconDense
        module procedure SolveQPDirectEqconSparse
        module procedure SolveQPDirectIneqconDense
        module procedure SolveQPDirectIneqconSparse
    end interface

    contains 

    !==================================================================!
    !                                                                  !
    !                              ROUTINES                            !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                             SOLVERS                              !
    !------------------------------------------------------------------!

    ! Classic unconstrained problem - easy gg no re
    subroutine SolveQPDirectNocon(hessJ, gradJ, x, flag)

        ! Description
        !============
        ! Solve an unconstrained quadratic problem of the form 
        !
        !   min_x 1/2 x' B x + g' x + c 
        !
        ! for the unknown vector x. The solution is obtained by 
        ! directly solving B x = -g, using intrinsic routines of the
        ! hessian approximation B rather than calling the solver 
        ! directly (this allows for more flexibility).

        ! Note: in the implementation, B is renamed to hessJ and g to
        ! gradJ. 
        
        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationUDT), intent(in)  :: hessJ 
        real(R8), intent(in)                        :: gradJ(:)
        real(R8), allocatable, intent(out)          :: x(:)
        integer(I8), intent(out)                    :: flag 

        ! Compute solution
        !=================
        x = hessJ%InverseHessianVectorProduct(-gradJ)

        ! Set flag to zero (might need to adjust later if product might
        ! not succeed, e.g. when inverting matrix etc)
        flag = 0

    end subroutine

    ! (linear) equality constrained problem, dense input
    subroutine SolveQPDirectEqconDense(hessJ, gradJ, jacG, b, x, lambda, flag)

        ! Description
        !============
        ! Solve a linear equality constrained quadratic problem of the
        ! form
        !
        !   min_x   0.5*x' hessJ x + gradJ' x 
        !   s.t.    jacG*x - b = 0
        !
        ! for the unknown variable x and the lagrange multipliers of
        ! the constraints lambda. It is assumed that jacG is in dense 
        ! format. The lagrangian is formulated as
        !
        !   L(x, lambda) = 0.5 x' hessJ x + gradJ' x + lambda' (jacG*x - b)
        !
        ! The solution follows from formulating the KKT conditions:
        !
        !   gradxL          = hessJ x + gradJ + lambda' jacG = 0
        !   gradlambdaL     = jacG*x - b = 0
        !
        ! Linearizing these conditions w.r.t. x and lambda yields
        !
        !   [hessJ, jacG'] [x]       = [-gradxL]
        !   [jacG,    0  ] [lambda]  = [b]
        !
        ! Solving this system yields the solution. Note that it is 
        ! linear in terms of x and lambda, and that for x0 = 0 and
        ! lambda0 = 0, gradxL reduces to gradJ (so we linearize
        ! around that point)

        ! The main issue is that we don't necessarily know hessJ 
        ! explicitly, so depending on the actual type of hessian 
        ! representation, we wil need a different algorithm:
        ! - sparse/dense representation: 'easy', just solve the 
        ! coupled system -> this is done here
        ! - limited-memory/others: hard -> this is not done here and
        ! an error will be thrown. Use a different algorithm (e.g. an 
        ! augmented lagrangian or other penalty method) to solve this. 

        ! Note: the parsed values for x and lambda serve as initial 
        ! guess (though that's unused here)

        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationUDT), intent(in)  :: hessJ 
        real(R8), intent(in)                        :: gradJ(:), &
            jacG(:, :), b(:)
        real(R8), intent(inout)                     :: x(size(gradJ)), lambda(size(b))

        ! Auxiliary
        integer(I8)                     :: flag, neq
        real(R8), allocatable           :: K(:, :), sol(:), rhs(:), &
            aux(:, :)

        ! Initialize
        !===========
        ! Get dimensions
        associate(nphi => hessJ%nphi) 
        neq = size(b)

        ! Check dimensions
        if (size(gradJ) /= nphi) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between cost function hessian and gradient')
        end if 
        if (size(jacG, 1) /= size(b)) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between constraint jacobian and constraint right hand side')
        end if 
        if (size(jacG, 2) /= nphi) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between constraint jacobian and cost function hessian')
        end if 

        ! Allocate 
        allocate(K(nphi+neq,nphi+neq), rhs(nphi+neq))

        ! Initialize
        K       = 0
        rhs(1:nphi) = -gradJ 
        rhs(nphi+1:nphi+neq) = b
        sol = rhs

        ! Compute
        !========
        ! Get dense hessian
        aux = hessJ%GetFullHessian()

        ! Construct the KKT system
        K(1:nphi, 1:nphi)           = aux
        K(1:nphi, nphi+1:nphi+neq)  = transpose(jacG)
        K(nphi+1:nphi+neq, 1:nphi)  = jacG 

        ! Solve the KKT system
        call SolveDenseLinearSystemDI(K, rhs, sol, flag)

        ! Check convergence
        if (flag /= 0) then 
            print *, 'SolveQPDirectEqcon: linear solver did not converge'
        end if 

        ! Unpack the solution
        x = sol(1:nphi)
        lambda = sol(nphi+1:nphi+neq)

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! (linear) equality constrained problem, sparse input
    subroutine SolveQPDirectEqconSparse(hessJ, gradJ, jacG, b, x, lambda, flag)

        ! Description
        !============
        ! Solve a linear equality constrained quadratic problem of the
        ! form
        !
        !   min_x   0.5*x' hessJ x + gradJ' x 
        !   s.t.    jacG*x - b = 0
        !
        ! for the unknown variable x and the lagrange multipliers of
        ! the constraints lambda. It is assumed that jacG is in sparse 
        ! format. The lagrangian is formulated as
        !
        !   L(x, lambda) = 0.5 x' hessJ x + gradJ' x + lambda' (jacG*x - b)
        !
        ! The solution follows from formulating the KKT conditions:
        !
        !   gradxL          = hessJ x + gradJ + lambda' jacG = 0
        !   gradlambdaL     = jacG*x - b = 0
        !
        ! Linearizing these conditions w.r.t. x and lambda yields
        !
        !   [hessJ, jacG'] [x]       = [-gradxL]
        !   [jacG,    0  ] [lambda]  = [b]
        !
        ! Solving this system yields the solution. Note that it is 
        ! linear in terms of x and lambda, and that for x0 = 0 and
        ! lambda0 = 0, gradxL reduces to gradJ (so we linearize
        ! around that point)

        ! The main issue is that we don't necessarily know hessJ 
        ! explicitly, so depending on the actual type of hessian 
        ! representation, we wil need a different algorithm:
        ! - sparse/dense representation: 'easy', just solve the 
        ! coupled system -> this is done here
        ! - limited-memory/others: hard -> this is not done here and
        ! an error will be thrown. Use a different algorithm (e.g. an 
        ! augmented lagrangian or other penalty method) to solve this. 

        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationUDT), intent(in)  :: hessJ 
        real(R8), intent(in)                        :: gradJ(:), &
            b(:)
        real(R8), intent(inout)                     :: x(size(gradJ)), lambda(size(b))
        type(MySparseUDT), intent(in)               :: jacG

        ! Auxiliary
        integer(I8)                     :: flag, neq
        real(R8), allocatable           :: sol(:), rhs(:)
        type(MySparseUDT)               :: K, aux

        ! Initialize
        !===========
        ! Get dimensions
        associate(nphi => hessJ%nphi) 
        neq = size(b)

        ! Check dimensions
        if (size(gradJ) /= nphi) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between cost function hessian and gradient')
        end if 
        if (jacG%nrow /= size(b)) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between constraint jacobian and constraint right hand side')
        end if 
        if (jacG%ncol /= nphi) then 
            call gdErrorHandler('SolveQPDirectEqcon: inconsistent dimensions ' // & 
                'between constraint jacobian and cost function hessian')
        end if 

        ! Allocate 
        allocate(rhs(nphi+neq))

        ! Initialize
        rhs(1:nphi) = -gradJ 
        rhs(nphi+1:nphi+neq) = b
        sol = rhs

        ! Compute
        !========
        ! Get sparse hessian
        aux = hessJ%GetSparseHessian()

        ! Construct the KKT system
        K = aux%Concatenate(jacG%Transpose(), 2)
        K = K%Concatenate(jacG%Concatenate(SpZeros(neq, neq), 2), 1)

        ! Solve the KKT system
        call SolveSparseLinearSystemDI(K, rhs, sol, flag)

        ! Check convergence
        if (flag /= 0) then 
            print *, 'SolveQPDirectEqcon: linear solver did not converge'
        end if 

        ! Unpack the solution
        x = sol(1:nphi)
        lambda = sol(nphi+1:nphi+neq)

        ! Housekeeping
        !=============
        end associate


    end subroutine

    ! (linear) inequality constrained problem, dense input
    subroutine SolveQPDirectIneqconDense(hessJ, gradJ, jacG, b, jacH, c, &
        x, lambda, mu, flag, maxit, tol, verbosity)

        ! Description
        !============
        ! Solver for the general linear inequality constrained problem 
        ! using a direct method. Note that the presence of inequality
        ! constraints immediately leads to a non-smooth problem, and 
        ! that we therefore need to iterate to determine the correct
        ! active set. We solve the general problem:
        !
        !       min_x 1/2   x' hessJ x + gradJ' x 
        !       s.t.        jacG x - b = 0
        !       s.t.        jacH x - c <= 0
        !   
        ! Defining l and m as the lagrange multipliers of the 
        ! constraints, the Lagrangian is formulated as
        !
        !       L = x' hessJ x + gradJ' x + l'(jacG x - b) + m_A'(jacH_A x - c_A)
        !
        ! where the subscript A denotes the indices of the active 
        ! constraints. In order to solve this system, we reformulate the
        ! inequality constraints using a nonlinear complementarity 
        ! function and iterate until the active set does not change
        ! anymore. For the ncp function, we take the max function and
        ! assume that alpha = 1

        ! Note: if only inequality constraints are present, pass 
        ! matrices jacG and b with correct 0(-by-nphi) dimensions. 

        ! Note: if the solver doesn't converge, this can be due to 
        ! insufficient number of iterations (flag equals -1) or due to
        ! linear solver issues during iteration (flag > 0). A message 
        ! will be shown if verbosity > 0

        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationUDT), intent(in)  :: hessJ 
        integer(I8), intent(in)                     :: maxit
        real(R8), intent(in)                        :: gradJ(:), &
            jacG(:, :), b(:), jacH(:, :), c(:), tol
        real(R8), intent(inout)                     :: x(size(gradJ)), &
            lambda(size(b)), mu(size(c))

        ! Auxiliary
        logical                         :: converged
        logical, allocatable            :: A(:)
        integer(I8)                     :: flag, verbosity
        real(R8), allocatable           :: K(:, :), sol(:), rhs(:), &
            aux(:, :), dncp(:, :), temp(:, :), ncp(:)

        ! Loop
        integer(I8)                     :: it, i

        ! Initialize
        !===========
        ! Check hessian representation
        aux = hessJ%GetFullHessian()

        ! Get problem dimensions
        associate(&
            nphi    => hessJ%nphi,      &
            neq     => size(jacG, 1),   &
            nineq   => size(jacH, 1))

        ! Initialize KKT matrix etc
        allocate(K(nphi+neq+nineq, nphi+neq+nineq), sol(nphi+neq+nineq), &
            dncp(nineq, nphi+neq+nineq))
        K = 0
        sol = 0
        rhs = sol

        ! Compute solution
        !=================
        ! Initialize
        converged = .false. 
        allocate(temp(nineq, 2))
        temp(:, 2) = 0

        ! Initialize constant parts of KKT matrix
        K = 0
        K(1:nphi, 1:nphi) = aux
        K(1:nphi, nphi+1:nphi+neq) = transpose(jacG)
        K(nphi+1:nphi+nineq, 1:nphi) = jacG 
        
        ! Print
        if (verbosity > 1) then 
            print *, 'Solving inequality constrained QP problem with ', maxit, ' iterations'
            print *, '| iterate | maxabs(dx, dlambda, dmu) | maxabs(res) | tol |'
        end if 

        ! Loop
        it = 1
        do while ((.not. converged) .and. (it <= maxit))

            ! Compute ncp function and active constraints
            temp(:, 1) = matmul(jacH, x) - c + mu
            ncp = maxval(temp, 2) - mu
            A = (temp(:, 1)) >= 0

            ! Update rhs
            rhs(1:nphi) = -(matmul(aux, x) + gradJ &
                + matmul(transpose(jacG), lambda) & 
                + matmul(transpose(jacH), mu))
            rhs(nphi+1:nphi+neq) = -(matmul(jacG, x) - b)
            rhs(nphi+neq+1:nphi+neq+nineq) = -(ncp)

            ! Print iterate 
            if (verbosity > 1) then 
                print *, '|', it,'|', maxval(abs(sol)),'|', maxval(abs(rhs)),'|', tol, '|' 
            end if
            
            ! Check convergence
            if (maxval(abs(rhs)) < tol) then 
                converged = .true.
                exit
            end if 

            ! Compute ncp linearization
            dncp = 0
            do i = 1, nineq 
                if (A(i)) then 
                    dncp(i, 1:nphi) = jacH(i, :)
                else
                    dncp(i, nphi+neq+i) = -1
                end if
            end do

            ! Update KKT matrix
            K(nphi+neq+1:nphi+neq+nineq, :) = dncp 
            K(1:nphi+neq, nphi+neq+1:nphi+neq+nineq) = transpose(dncp(:, 1:nphi+neq))

            ! Compute update
            call SolveDenseLinearSystemDI(K, rhs, sol, flag)

            ! Check the flag
            if (flag /= 0) then 
                ! Quit the loop
                converged = .true. 
                exit
            end if 

            ! Assign update
            x = x + sol(1:nphi)
            lambda = lambda + sol(nphi+1:nphi+neq)
            mu = mu + sol(nphi+neq+1:nphi+neq+nineq)

            ! Project mu
            where(mu < 0) mu = 0
            
            ! Update iteration counter
            it = it + 1

        end do

        ! Checks
        !=======
        ! If not converged, check why
        if (verbosity > 0) then 
            if (converged .and. (flag /= 0)) then 
                ! Solver did not converge at some point
                print *, 'SolveQPDirectIneqconDense: linear solver did not converge, exited prematurely'
            end if 
            if (.not. converged) then 
                ! Maximum number of iterations reached
                print *, 'SolveQPDirectIneqconDense: linear solver did not converge, exited prematurely'
            end if 
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

    ! (linear) inequality constrained problem, sparse input
    subroutine SolveQPDirectIneqconSparse(hessJ, gradJ, jacG, b, jacH, c, &
        x, lambda, mu, flag, maxit, tol, verbosity)

        ! Description
        !============
        ! Solver for the general linear inequality constrained problem 
        ! using a direct method. Note that the presence of inequality
        ! constraints immediately leads to a non-smooth problem, and 
        ! that we therefore need to iterate to determine the correct
        ! active set. We solve the general problem:
        !
        !       min_x 1/2   x' hessJ x + gradJ' x 
        !       s.t.        jacG x - b = 0
        !       s.t.        jacH x - c <= 0
        !   
        ! Defining l and m as the lagrange multipliers of the 
        ! constraints, the Lagrangian is formulated as
        !
        !       L = x' hessJ x + gradJ' x + l'(jacG x - b) + m_A'(jacH_A x - c_A)
        !
        ! where the subscript A denotes the indices of the active 
        ! constraints. In order to solve this system, we reformulate the
        ! inequality constraints using a nonlinear complementarity 
        ! function and iterate until the active set does not change
        ! anymore. For the ncp function, we take the max function and
        ! assume that alpha = 1

        ! Note: if only inequality constraints are present, pass 
        ! matrices jacG and b with correct 0(-by-nphi) dimensions. 

        ! Note: if the solver doesn't converge, this can be due to 
        ! insufficient number of iterations (flag equals -1) or due to
        ! linear solver issues during iteration (flag > 0). A message 
        ! will be shown if verbosity > 0

        ! Declare variables
        !==================
        ! Arguments
        class(HessianApproximationUDT), intent(in)  :: hessJ 
        integer(I8), intent(in)                     :: maxit
        real(R8), intent(in)                        :: gradJ(:), &
            b(:), c(:), tol
        real(R8), intent(inout)                     :: x(size(gradJ)), &
            lambda(size(b)), mu(size(c))
        type(MySparseUDT), intent(in)               :: jacG, jacH

        ! Auxiliary
        logical                         :: converged
        logical, allocatable            :: A(:)

        integer(I8)                     :: flag, verbosity, nic
        integer(I8), allocatable        :: conind(:), rowI(:), colI(:) 

        real(R8), allocatable           :: sol(:), rhs(:), &
            temp(:, :), ncp(:), valI(:)

        type(MySparseUDT)               :: aux, dncp, K, Kc, dncpI, &
            jacGT, jacHT

        ! Loop
        integer(I8)                     :: it, j

        ! Initialize
        !===========
        ! Check hessian representation
        aux = hessJ%GetSparseHessian()

        ! Get problem dimensions
        associate(&
            nphi    => hessJ%nphi,      &
            neq     => jacG%nrow,       &
            nineq   => jacH%nrow)

        ! Initialize 
        allocate(sol(nphi+neq+nineq))
        sol = 0
        rhs = sol
        conind = [(j, j = 1, nineq)]
        jacGT = jacG%Transpose()
        jacHT = jacH%Transpose()

        ! Compute solution
        !=================
        ! Initialize
        converged = .false. 
        allocate(temp(nineq, 2))
        temp(:, 2) = 0

        ! Initialize constant parts of KKT matrix
        Kc = aux%Concatenate(jacG%Transpose(), 2)
        Kc = Kc%Concatenate(jacG%Concatenate(SpZeros(neq, neq), 2), 1)
        
        ! Print
        if (verbosity > 1) then 
            print *, 'Solving inequality constrained QP problem with ', maxit, ' iterations'
            print *, '| iterate | maxabs(dx, dlambda, dmu) | maxabs(res) | tol |'
        end if 

        ! Loop
        it = 1
        do while ((.not. converged) .and. (it <= maxit))

            ! Compute ncp function and active constraints
            temp(:, 1) = jacH%MatrixVectorProduct(x) - c + mu
            ncp = maxval(temp, 2) - mu
            A = (temp(:, 1)) >= 0

            ! Update rhs
            rhs(1:nphi) = -(aux%MatrixVectorProduct(x) + gradJ &
                + jacGT%MatrixVectorProduct(lambda) & 
                + jacHT%MatrixVectorProduct(mu))
            rhs(nphi+1:nphi+neq) = -(jacG%MatrixVectorProduct(x) - b)
            rhs(nphi+neq+1:nphi+neq+nineq) = -(ncp)

            ! Print iterate 
            if (verbosity > 1) then 
                print *, '|', it,'|', maxval(abs(sol)),'|', maxval(abs(rhs)),'|', tol, '|' 
            end if

            ! Check convergence
            if (maxval(abs(rhs)) < tol) then 
                converged = .true.
                exit
            end if 

            ! Compute ncp linearization
            nic = count(.not. A) ! number of inactive constraints
            dncp = jacH%SetZeroRows(.not. A) ! set inactive to zero for now
            dncp = dncp%Concatenate(SpZeros(nineq, neq), 2) ! add linearization w.r.t. lambda - zero
            allocate(rowI(nic), colI(nic), valI(nic))
            rowI = pack(conind, .not. A)
            colI = pack(conind, .not. A) 
            valI = -1.0
            dncpI = ConstructMySparse(rowI, colI, valI, nineq, nineq)
            deallocate(rowI, colI, valI)

            ! Update KKT matrix
            K = Kc%Concatenate(dncp%Transpose(), 2)
            K = K%Concatenate(dncp%Concatenate(dncpI, 2), 1)

            ! Compute update
            call SolveSparseLinearSystemDI(K, rhs, sol, flag)

            ! Check the flag
            if (flag /= 0) then 
                ! Quit the loop
                converged = .true. 
                exit
            end if 

            ! Assign update
            x = x + sol(1:nphi)
            lambda = lambda + sol(nphi+1:nphi+neq)
            mu = mu + sol(nphi+neq+1:nphi+neq+nineq)

            ! Project mu
            where(mu < 0) mu = 0

            ! Update iteration counter
            it = it + 1

        end do

        ! Checks
        !=======
        ! If not converged, check why
        if (verbosity > 0) then 
            if (converged .and. (flag /= 0)) then 
                ! Solver did not converge at some point
                print *, 'SolveQPDirectIneqconSparse: linear solver did not converge, exited prematurely'
            end if 
            if (.not. converged) then 
                ! Maximum number of iterations reached
                print *, 'SolveQPDirectIneqconSparse: linear solver did not converge, exited prematurely'
            end if 
        end if 

        ! Housekeeping
        !=============
        end associate

    end subroutine

end module