!==================================================================!
!                                                                  !
!                        DOCUMENTATION                             !
!                                                                  !
!==================================================================!
! This module provides the graph type and its methods. Currently, this 
! contains only a crude homegrown implementation of the most important
! graph functionality for grid generation (e.g. determination of 
! number of subgraphs, basic interconnection checking etc). In the 
! future, it is likely better to couple to a dedicated optimized 
! graph handling package, which is why all implementation here is 
! private, including the type's data (this way, code that builds 
! upon this module shouldn't be impacted by changes under the hood). 
! Some data we have public as we know this data won't change (e.g. 
! number of nodes, connections per node, etc.)

! We support a number of graph types, depending on desired properties
! that can be exploited. Therefore, the main type is the most general,
! directed, multiply connected (duplicate edges), looped (edges of a
! single vertex) graph type there is. Consequently, there is not much
! exploitation of structure available, nor can some queries not be 
! performed. To this end, it is typicaly better to use one of the 
! inheriting types that have more efficient methods (and sometimes 
! additional functionality). Do note that in order to use additional 
! functionality, one has to specifically declare the type as a class of
! at least the type desired (or use the 'select type' construct)

! IMPORTANT: in the general graph implementation, the vertex IDs are
! mapped to a 1-nv array (edges are also adjusted accordingly) by using 
! a binary search method for sorted arrays. This may not be very 
! performant for large graphs and may at some point be replaced by 
! a hash map (if such implementation ever becomes available). 

module mod_graph
    ! Load modules
    use mod_precision
    use mod_errorhandler
    use mod_sort
    use mod_dynamicarrays

    implicit none
    private ! default

    !==================================================================!
    !                                                                  !
    !                            TYPES                                 !
    !                                                                  !
    !==================================================================!

    ! Main graph type
    type, public :: GraphUDT

        ! Description
        !============
        ! Main graph type, most general type description possible. The 
        ! graph simply consists of a set of vertices (v) and edges (e),
        ! where the latter are assumed to be directed (i.e. an edge is 
        ! actually an arc, but ok). Edges may have multiplicity etc, 
        ! vertices cannot. (so this is actually a multidigraph)

        ! We define some basic derived quantities that should be useful
        ! for any graph (e.g. the edges of a vertex). We assume that 
        ! these quantities should only be updated when vertices/edges 
        ! are deleted.

        ! Note: the edges do not store the actual vertex IDs, but the 
        ! index into v to get the correct vertex ID. It is assumed that 
        ! the vertex index always goes from 1 to graph%nv, so when 
        ! manipulating graphs, care must be taken to update edges etc
        ! accordingly. 

        ! Note: the vertices will be sorted according to their original
        ! vertex ID (as a side effect of the unique function). This 
        ! will be exploited when doing mapping operations. 

        ! Define basic quantities and sizes
        integer(I8)                             :: nv, ne 
        integer(I8), allocatable, dimension(:)  :: v, ev1, ev2

        ! Define derived quantities 
        integer, allocatable, dimension(:)      :: ve, vep1, vep2 

    contains

        ! Basics
        !=======
        ! Constructor 
        procedure :: Construct          => ConstructGraph 

        ! Derived quantity computation
        procedure :: UpdateInterconnectionData      => UpdateGraphInterconnectionData

        ! Edge deletion
        procedure, private :: DeleteEdgeVertexBased     => DeleteGraphEdgeVertexBased
        procedure, private :: DeleteEdgeLogicalBased    => DeleteGraphEdgeLogicalBased
        generic :: DeleteEdges  => DeleteEdgeVertexBased, DeleteEdgeLogicalBased

        ! Vertex deletion
        procedure, private :: DeleteGraphVertexIndexBased
        procedure, private :: DeleteGraphVertexLogicalBased
        generic :: DeleteVertex => DeleteGraphVertexIndexBased, &
            DeleteGraphVertexLogicalBased

        ! Getters
        procedure :: GetVertexIndex     
        procedure :: GetVertexEdges     => GetGraphVertexEdges
        procedure :: GetVertexEdgesDirectional  => GetGraphVertexEdgesDirectional

        ! Operators
        !==========
        ! Flood
        procedure :: Flood      => FloodGraph 

        ! Condensation
        procedure :: Condense   => CondenseGraph

        ! Cycle checking

    end type

    type, extends(GraphUDT), public :: UGraphUDT 
        
        ! Description
        !============
        ! This graph type represents a general yet undirected graph. 
        ! Most functionality is taken over from the parent type, though
        ! some procedures are slightly modified, since we do not store
        ! each edge double (i.e. the edge [v1, v2] represents the 
        ! connection v1, v2 and v2, v1)

        ! Note: the directional edge getters are still supported, but 
        ! should be used wisely

    contains

        ! Basics
        !=======
        ! Constructor (unchanged)

        ! Derived quantity computation (unchanged)

        ! Edge deletion 
        procedure, private :: DeleteEdgeVertexBased     => DeleteUGraphEdgeVertexBased

        ! Getters (unchanged)

        ! Operators
        !==========
        ! Flood
        procedure :: Flood      => FloodUGraph 

        ! Condensation
        procedure :: Condense   => CondenseUGraph

        ! Connection 
        procedure :: IsConnected    => IsUGraphConnected

    end type

contains     

    !==================================================================!
    !                                                                  !
    !                           ROUTINES                               !
    !                                                                  !
    !==================================================================!

    !------------------------------------------------------------------!
    !                         GENERAL GRAPH                            !
    !------------------------------------------------------------------!
    ! Constructor
    subroutine ConstructGraph(graph, ev1, ev2, v)

        ! Description
        !============
        ! Constructor for the general graph. No assumptions on the 
        ! underlying structure are made, and edges are assumed to be 
        ! arcs (i.e. there's a sense of direction). Vertices are made 
        ! unique and sorted.

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                         :: graph 
        integer(I8), dimension(:), intent(in)   :: ev1, ev2, v

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: vu

        ! Loop
        integer(I8)                             :: i, ind 

        ! Checks
        !=======
        if (size(ev1) /= size(ev2)) then 
            call gdErrorHandler('ConstructGraph: inconsistent size of ' // &
                'input arguments')
        end if 

        ! Initialize
        !===========
        ! Compute vertices
        call Unique(v, vu) 
        graph%v = vu 
        graph%nv = size(vu)

        ! Compute edges
        graph%ne = size(ev1)
        graph%ev1 = ev1 
        graph%ev2 = ev2
        do i = 1, size(ev1)
            ! First edge vertex
            ind = graph%GetVertexIndex(ev1(i))
            if (ind == 0) then 
                call gdErrorHandler('ConstructGraph: edge vertex is ' // &
                    'not given in vertex set, check input')
            end if 
            graph%ev1(i) = ind 

            ! Second edge vertex
            ind = graph%GetVertexIndex(ev2(i))
            if (ind == 0) then 
                call gdErrorHandler('ConstructGraph: edge vertex is ' // &
                    'not given in vertex set, check input')
            end if 
            graph%ev2(i) = ind 
        end do 
        

    end subroutine

    ! Interconnection data computation
    subroutine UpdateGraphInterconnectionData(graph)

        ! Description 
        !============
        ! Update the derived data that is based on the basic vertex and
        ! edge data. Note that interconnection data is recomputed from
        ! scratch, which may not be the most efficient approach 
        ! depending on the use case. 

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 

        ! Auxiliary 

        ! Compute interconnections
        !=========================
        ! Vertex edges


    end subroutine

    ! Edge deletion (vertex based)
    subroutine DeleteGraphEdgeVertexBased(graph, v1, v2)

        ! Description
        !============
        ! Delete the edge(s) with vertex indices [v1, v2]. 
        ! If no edges were found, no error is thrown and the 
        ! graph remains unchanged. Edges with vertices [v2, v1] will not
        ! be deleted, as it is assumed the graph is directional. 

        ! Note: the vertex indices here are local indices, so not equal
        ! to those in v!

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 
        integer(I8), intent(in)             :: v1, v2

        ! Auxiliary
        logical, allocatable, dimension(:)  :: delvec 

        ! Initialize
        !===========
        ! Allocate
        allocate(delvec(graph%ne))
        delvec = .false. 

        ! Delete
        !=======
        ! Check edges with [v1, v2]
        where (v1 == graph%ev1 .and. v2 == graph%ev2) delvec = .true. 

        ! Delete edges (but not vertices)
        call graph%DeleteEdgeLogicalBased(delvec)
        
    end subroutine

    ! Edge deletion (logical based)
    subroutine DeleteGraphEdgeLogicalBased(graph, delvec)

        ! Description
        !============
        ! Delete edges that are 'true' in delvec

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 
        logical, dimension(:), intent(in)   :: delvec 

        ! Check
        !======
        if (size(delvec) /= graph%ne) then 
            call gdErrorHandler('DeleteGraphEdgeLogicalBased: ' // & 
                'deletion vector has incompatible dimension')
        end if

        ! Delete
        !=======
        graph%ev1 = pack(graph%ev1, .not. delvec)
        graph%ev2 = pack(graph%ev2, .not. delvec)
        graph%ne = size(graph%ev1)

    end subroutine

    ! Vertex deletion (index based)
    subroutine DeleteGraphVertexIndexBased(graph, vertind)

        ! Description
        !============
        ! Remove a single vertex from the graph

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 
        integer(I8), intent(in)             :: vertind

        ! Auxiliary
        integer(I8), allocatable, dimension(:)  :: eind
        logical, allocatable, dimension(:)      :: deledge
        
        ! Loop
        integer(I8)                             :: i 


        ! Check
        !======
        if (vertind > graph%nv .or. vertind < 0) then 
            call gdErrorHandler('DeleteGraphVertexIndexBased: ' // & 
                'deletion index is out of bounds')
        end if

        ! Mark for deletion
        !==================
        ! Initialize
        allocate(deledge(graph%ne))
        deledge = .false. 

        ! Get vertex edge indices
        eind = graph%GetVertexEdges(i)

        ! Mark for deletion
        deledge(eind) = .true. 

        ! Delete
        !=======
        graph%v = [graph%v(1:vertind-1), graph%v(vertind+1:)]
        graph%ev1 = pack(graph%ev1, .not. deledge)
        graph%ev2 = pack(graph%ev2, .not. deledge)
        graph%nv = size(graph%v)
        graph%ne = size(graph%ev1)
        where (graph%ev1 > vertind) graph%ev1 = graph%ev1 - 1
        where (graph%ev2 > vertind) graph%ev2 = graph%ev2 - 1

    end subroutine 

    ! Vertex deletion (logical based)
    subroutine DeleteGraphVertexLogicalBased(graph, delvec)

        ! Description
        !============
        ! Delete a graph vertex based on the logical 'delvec'. Any edges
        ! connecting to this vertex are also deleted. 

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 
        logical, dimension(:), intent(in)   :: delvec

        ! Auxiliary
        integer(I8)                             :: vind 
        integer(I8), allocatable, dimension(:)  :: eind, mapv 
        logical, allocatable, dimension(:)      :: deledge
        
        ! Loop
        integer(I8)                             :: i 

        ! Check
        !======
        if (size(delvec) /= graph%nv) then 
            call gdErrorHandler('DeleteGraphVertexLogicalBased: ' // & 
                'deletion vector has incompatible dimension')
        end if

        ! Mark for deletion
        !==================
        allocate(deledge(graph%ne))
        deledge = .false. 
        do i = 1, graph%nv 
            if (delvec(i)) then 
                ! Get vertex edge indices
                eind = graph%GetVertexEdges(i)

                ! Mark for deletion
                deledge(eind) = .true. 
            end if 
        end do 

        ! Construct mapping
        allocate(mapv(graph%nv))
        mapv = 0
        vind = 0
        do i = 1, graph%nv
            if (delvec(i)) then 
                cycle 
            end if 
            vind = vind + 1
            mapv(i) = vind 
        end do 

        ! Delete
        !=======
        graph%v = pack(graph%v, .not. delvec)
        graph%ev1 = pack(graph%ev1, .not. deledge)
        graph%ev2 = pack(graph%ev2, .not. deledge)
        graph%nv = size(graph%v)
        graph%ne = size(graph%ev1)
        graph%ev1 = mapv(graph%ev1)
        graph%ev2 = mapv(graph%ev2)

    end subroutine

    ! Vertex index getter
    function GetVertexIndex(graph, vID) result(ind)

        ! Description
        !============
        ! This function returns the local vertex index of a vertex 
        ! with ID 'vID' in the v array of the graph. If the vertex 
        ! is not present, zero is returned. 

        ! Note: this is just a wrapper for any search algorithm to 
        ! provide proper encapsulation. 

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                 :: graph 
        integer(I8), intent(in)         :: vID 
        integer(I8)                     :: ind 

        ! Search
        !=======
        ind =  SearchSortedArray(graph%v, vID)

    end function

    ! Vertex edge getter
    function GetGraphVertexEdges(graph, ind) result(eind)

        ! Description
        !============
        ! Get the edges of a vertex in index format 
        ! Note: since the basic graph object does not contain any 
        ! interconnection data, we need to compute this from scratch. 
        ! This may not be the best implementation and is probably best
        ! overwritten by inheriting types with better data 
        ! structures (e.g. that keep a vertex to edge pointer)

        ! Note: these are all the vertex's edges, including ingoing
        ! and outgoing! 

        ! Declare variables
        !==================
        ! Declare variables
        class(GraphUDT)                         :: graph 
        integer(I8), intent(in)                 :: ind 
        integer(I8), allocatable, dimension(:)  :: eind

        ! Auxiliary
        logical, allocatable, dimension(:)      :: hasvert

        ! Loop
        integer(I8)                             :: k 

        ! Find edges
        !===========
        hasvert = graph%ev1 == ind .or. graph%ev2 == ind 
        allocate(eind(count(hasvert)))
        eind = pack([(k, k = 1, graph%ne)], hasvert)

    end function

    ! Vertex edge getter, directional
    function GetGraphVertexEdgesDirectional(graph, ind, outwards) &
        result(eind)

        ! Description
        !============
        ! This function gets the edges of a vertex that either go 
        ! outwards (so edges with ind in ev1) or inwards (edges with 
        ! ind in ev2). This is checked with the 'outwards' logical. 

        ! Declare variables
        !==================
        ! Declare variables
        class(GraphUDT)                         :: graph 
        integer(I8), intent(in)                 :: ind 
        integer(I8), allocatable, dimension(:)  :: eind
        logical, intent(in)                     :: outwards 

        ! Auxiliary
        logical, allocatable, dimension(:)      :: hasvert

        ! Loop
        integer(I8)                             :: k 

        ! Find edges
        !===========
        if (outwards) then 
            hasvert = graph%ev1 == ind 
        else
            hasvert = graph%ev2 == ind 
        end if  
        eind = pack([(k, k = 1, graph%ne)], hasvert)
        
    end function

    ! Flood
    function FloodGraph(graph, startvertind) result(vertind)

        ! Description
        !============
        ! This function 'floods' the graph, starting from the vertex 
        ! index given in 'startvertind'. A 'flood' means that all vertices 
        ! that are somehow connected by edges to this vertex are found.
        ! If the graph is fully connected, all entries of vertind will
        ! be true (this indicates if a vertex is present in the flooded
        ! part of the graph). Otherwise, the graph exists of multiple 
        ! non-connected parts, called 'subgraphs', which each can be 
        ! found by flooding from a non-flooded vertex. 

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT)                     :: graph 
        integer(I8), intent(in)             :: startvertind
        logical, allocatable, dimension(:)  :: vertind 

        ! Auxiliary
        integer(I8)                             :: nextv
        integer(I8), allocatable, dimension(:)  :: ve
        logical, allocatable, dimension(:)      :: istraceable 

        ! Initialize
        !===========
        ! Check
        if ((startvertind > graph%nv) .or. (startvertind < 1)) then 
            call gdErrorHandler('FloodGraph: start vertex index is out ' // & 
                'of bounds')
        end if 

        ! Initialize
        allocate(vertind(graph%nv),  istraceable(graph%nv))
        vertind = .false. 
        istraceable = .false. 

        ! Flood
        !======
        ! Better call Master Chief
        nextv = startvertind
        vertind(nextv) = .true. 
        do while (nextv /= 0)
            ! Get the vertex outgoing edges
            ve = graph%GetVertexEdgesDirectional(nextv, .true.)

            ! Set vertices as traceable
            istraceable([graph%ev1(ve), graph%ev2(ve)]) = .true. 

            ! Set vertices that were already found as non-traceable
            where (vertind) istraceable = .false.  

            ! Get the next vertex index 
            nextv = findloc(istraceable, .true., 1, back=.false.)

        end do 

    end function 

    ! Graph condensation
    subroutine CondenseGraph(graph)

        ! Description
        !============
        ! This routine condenses a graph by removing any duplicate 
        ! edges and any edges that consist of the same vertex IDs. 
        ! Dangling/pendant vertices are not removed
        ! of course. 

        ! Declare variables
        !==================
        ! Arguments
        class(GraphUDT), intent(inout)              :: graph 

        ! Auxiliary
        integer(I8), allocatable, dimension(:)      :: ev1, ev2, &
            sortind, dev1, dev2 
        logical, allocatable, dimension(:)          :: delvec 

        ! Loop
        integer(I8)                                 :: si, ei

        ! Initialize
        !===========
        ! Allocate
        allocate(delvec(size(ev1)))
        delvec = .false. 

        ! Copy
        ev1 = graph%ev1
        ev2 = graph%ev2

        ! Sort edges rowwise if it is a non-directional graph to 
        ! capture similar but 'reversed' edges
        !if (.not. isdirectional) then 
        !    ! Sort
        !    do i = 1, size(ev1)
        !        if (ev1(i) > ev2(i)) then 
        !            ! swap 
        !            tmp = ev2(i)
        !            ev2(i) = ev1(i)
        !            ev1(i) = tmp 
        !        end if 
        !    end do 
        !end if 

        ! Sort edges 
        !===========
        ! Sort along ev1 
        allocate(sortind(size(ev1))) 
        call Sort(ev1, ind=sortind, ascend=.true.)
        ev2 = ev2(sortind)
        deallocate(sortind)

        ! Sort equal values of ev1 for ev2
        si = 0
        ei = 0
        do while (.true.)
            ! Get start index
            si = ei + 1

            ! Get end index
            ei = si 
            do while (ev1(si) == ev1(ei))
                ei = ei + 1
                if (ei > size(ev1)) then 
                    exit 
                end if 
            end do 

            ! Adjust due to loop
            ei = ei - 1 

            ! Sort this segment according to ev2
            call Sort(ev2(si:ei))
            
        end do 

        ! Remove duplicate edges and loops
        !=======================
        ! Compute differences
        dev1 = ev1(2:) - ev1(1:size(ev1))
        dev2 = ev2(2:) - ev2(1:size(ev2))

        ! Update deletion vector
        delvec = delvec .or. ([.false., dev1 == 0 .and. dev2 == 0])

        ! Remove loops
        delvec = delvec .or. (ev1 == ev2)

        ! Reconstruct graph
        !==================
        ! Reconstruct graph edges
        graph%ne = count(.not. delvec)
        graph%ev1 = pack(ev1, .not. delvec) 
        graph%ev2 = pack(ev2, .not. delvec)


    end subroutine

    !------------------------------------------------------------------!
    !                    GENERAL UNDIRECTED GRAPH                      !
    !------------------------------------------------------------------!

    ! Edge deletion (vertex based)
    subroutine DeleteUGraphEdgeVertexBased(graph, v1, v2)

        ! Description
        !============
        ! Delete the edge(s) with vertex indices [v1, v2]. 
        ! If no edges were found, no error is thrown and the 
        ! graph remains unchanged. Edges with vertices [v2, v1] will also
        ! be deleted, as it is assumed the graph is undirectional. 

        ! Note: the vertex indices here are local indices, so not equal
        ! to those in v!

        ! Declare variables
        !==================
        ! Arguments
        class(UGraphUDT)                    :: graph 
        integer(I8), intent(in)             :: v1, v2

        ! Auxiliary
        logical, allocatable, dimension(:)  :: delvec 

        ! Initialize
        !===========
        ! Allocate
        allocate(delvec(graph%ne))
        delvec = .false. 

        ! Delete
        !=======
        ! Check edges with [v1, v2]
        where (v1 == graph%ev1 .and. v2 == graph%ev2) delvec = .true.

        ! Check edges with [v2, v1]
        where (v2 == graph%ev1 .and. v1 == graph%ev2) delvec = .true. 

        ! Delete edges (but not vertices)
        call graph%DeleteEdgeLogicalBased(delvec)
        
    end subroutine

    ! Flood
    function FloodUGraph(graph, startvertind) result(vertind)

        ! Description
        !============
        ! This function 'floods' the graph, starting from the vertex 
        ! index given in 'startvertind'. A 'flood' means that all vertices 
        ! that are somehow connected by edges to this vertex are found.
        ! If the graph is fully connected, all entries of vertind will
        ! be true (this indicates if a vertex is present in the flooded
        ! part of the graph). Otherwise, the graph exists of multiple 
        ! non-connected parts, called 'subgraphs', which each can be 
        ! found by flooding from a non-flooded vertex. 

        ! Declare variables
        !==================
        ! Arguments
        class(UGraphUDT)                    :: graph 
        integer(I8), intent(in)             :: startvertind
        logical, allocatable, dimension(:)  :: vertind 

        ! Auxiliary
        integer(I8)                             :: nextv
        integer(I8), allocatable, dimension(:)  :: ve
        logical, allocatable, dimension(:)      :: istraceable 

        ! Initialize
        !===========
        ! Check
        if ((startvertind > graph%nv) .or. (startvertind < 1)) then 
            call gdErrorHandler('FloodGraph: start vertex index is out ' // & 
                'of bounds')
        end if 

        ! Initialize
        allocate(vertind(graph%nv),  istraceable(graph%nv))
        vertind = .false. 
        istraceable = .false. 

        ! Flood
        !======
        ! Better call Master Chief
        nextv = startvertind
        vertind(nextv) = .true. 
        do while (nextv /= 0)
            ! Get the vertex edges
            ve = graph%GetVertexEdges(nextv)

            ! Set vertices as traceable
            istraceable([graph%ev1(ve), graph%ev2(ve)]) = .true. 

            ! Set vertices that were already found as non-traceable
            where (vertind) istraceable = .false.  

            ! Get the next vertex index 
            nextv = findloc(istraceable, .true., 1, back=.false.)

        end do 

    end function 

    ! Graph condensation
    subroutine CondenseUGraph(graph)

        ! Description
        !============
        ! This routine condenses a graph by removing any duplicate 
        ! edges and any edges that consist of the same vertex IDs. 
        ! Dangling/pendant vertices are not removed
        ! of course. 

        ! Declare variables
        !==================
        ! Arguments
        class(UGraphUDT), intent(inout)             :: graph 

        ! Auxiliary
        integer(I8)                                 :: tmp 
        integer(I8), allocatable, dimension(:)      :: ev1, ev2, &
            sortind, dev1, dev2 
        logical, allocatable, dimension(:)          :: delvec 

        ! Loop
        integer(I8)                                 :: i, si, ei

        ! Initialize
        !===========
        ! Allocate
        allocate(delvec(size(ev1)))
        delvec = .false. 

        ! Copy
        ev1 = graph%ev1
        ev2 = graph%ev2

        ! Sort edges rowwise 
        do i = 1, size(ev1)
            if (ev1(i) > ev2(i)) then 
                ! swap 
                tmp = ev2(i)
                ev2(i) = ev1(i)
                ev1(i) = tmp 
            end if 
        end do 

        ! Sort edges 
        !===========
        ! Sort along ev1 
        allocate(sortind(size(ev1))) 
        call Sort(ev1, ind=sortind, ascend=.true.)
        ev2 = ev2(sortind)
        deallocate(sortind)

        ! Sort equal values of ev1 for ev2
        si = 0
        ei = 0
        do while (.true.)
            ! Get start index
            si = ei + 1

            ! Get end index
            ei = si 
            do while (ev1(si) == ev1(ei))
                ei = ei + 1
                if (ei > size(ev1)) then 
                    exit 
                end if 
            end do 

            ! Adjust due to loop
            ei = ei - 1 

            ! Sort this segment according to ev2
            call Sort(ev2(si:ei))
            
        end do 

        ! Remove duplicate edges and loops
        !=======================
        ! Compute differences
        dev1 = ev1(2:) - ev1(1:size(ev1))
        dev2 = ev2(2:) - ev2(1:size(ev2))

        ! Update deletion vector
        delvec = delvec .or. ([.false., dev1 == 0 .and. dev2 == 0])

        ! Remove loops
        delvec = delvec .or. (ev1 == ev2)

        ! Reconstruct graph
        !==================
        ! Reconstruct graph edges
        graph%ne = count(.not. delvec)
        graph%ev1 = pack(ev1, .not. delvec) 
        graph%ev2 = pack(ev2, .not. delvec)


    end subroutine

    ! Graph connectedness check
    function IsUGraphConnected(graph) result(isconnected)

        ! Description
        !============
        ! This function checks whether the graph is connected by 
        ! doing a flood operation and checking if each vertex was 
        ! found.

        ! Declare variables
        !==================
        ! Arguments
        class(UGraphUDT)                    :: graph 
        logical                             :: isconnected 

        ! Auxiliary
        logical, allocatable, dimension(:)  :: vertind 

        ! Compute
        !========
        ! Initialize
        isconnected = .false. 

        ! Hedge for trivial case
        if (graph%nv == 0) then 
            isconnected = .true.
            return  
        end if 

        ! Compute
        vertind = graph%Flood(1)
        if (all(vertind)) then 
            isconnected = .true. 
        end if 

    end function

end module
