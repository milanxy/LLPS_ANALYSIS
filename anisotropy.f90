!*************************************************************************************!

!Precision_control 
!*************************************************************************************!    
      module Precision_Control
      implicit none
      integer, parameter :: single = selected_real_kind(p=6,r=37) !single precision
      integer, parameter :: double = selected_real_kind(p=15,r=37) !double precision
      integer, parameter :: quad= selected_real_kind(p=33,r=4931)  !quadruple precision
      save
      end module


!*************************************************************************************!

!Parameter definition
!*************************************************************************************! 

      module Param_Control
      Use Precision_Control
      implicit none
      integer:: step_num, Nstep,Nstepskip, Npolymer, Nresi, Nbead_per_resi, Ndim
      real(kind=double) ::rho

      CHARACTER(LEN=100),PARAMETER ::FMT1="(3f8.3)"
      save
      end module

!*************************************************************************************!

!Variable definition
!*************************************************************************************! 
      module Var_Control
      Use Precision_Control
      Use Param_Control
      implicit none
      integer ::istep
      real(kind=double), allocatable, dimension(:) :: x, y, z
      real(kind=double)::boxlen, binsize
      real(kind=double), allocatable, dimension(:,:) :: com_lys,com_glu,com
      real(kind=double), allocatable, dimension(:,:) ::com_x,com_y,com_z
      real(kind=double), allocatable, dimension(:,:) ::com_droplet_x,com_droplet_y, com_droplet_z
      integer::Nbin
      real(kind=double),allocatable,dimension(:):: unnorm_summed_numbers


      integer, allocatable,dimension(:) ::  num_neighbor
      integer, allocatable,dimension(:,:)::neighborlist
      integer, allocatable,dimension(:) ::cluster_head, cluster_tail, cluster_size, cluster_list, cluster_checked
      integer,allocatable, dimension(:) :: state
      integer ::clustersum, av_droplet
      integer,allocatable,dimension(:,:)::listed_id,id_droplet,id_free_single
      integer,allocatable,dimension(:,:) ::id_particle_in_cluster
      real(kind=double):: div,av_div,av_sq_div
      save 
      end module

!*************************************************************************************!

!CORE
!*************************************************************************************! 
      program main
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::  iresi,ibin,ip
      integer ::ijunk
      real(kind=double) :: constant, vol, pi =3.14, num_ideal
      call INITIALIZE
      
      read(1,*)
      read(1,*)
      do iresi = 1 , Nresi
         read(1,*)
      end do
      
      do ijunk = 1, Nstepskip*4002
         read(1,*)
      end do
      
      do istep =1 , (Nstep-Nstepskip)
         
         call READ_FILE
         call FINDCOM
         call RDF
         call GETCLUSTER
         call GETANISOTROPY
         av_div = av_div + div
         av_sq_div=av_sq_div+div**2
      end do
      
      constant = (1.33333333333)* pi
      do ibin = 1, Nbin
         vol = constant*(((ibin*binsize)**3)-((ibin-1)*binsize)**3)
         num_ideal = vol*rho
         
         unnorm_summed_numbers(ibin)=unnorm_summed_numbers(ibin)/(num_ideal)
      end do
      do ibin = 1, Nbin
         !write(200,*) unnorm_summed_numbers(ibin)
         !write(100, *)ibin*binsize,unnorm_summed_numbers(ibin)/(Npolymer*(Nstep-Nstepskip))
      end do
      !write(1500,*) real(av_droplet)/(Nstep-Nstepskip)

      !do ip =1, Npolymer
      !   write(2000,*) id_droplet(ip,Nstep-Nstepskip)
      !end do
      !call GETDIFFUSION
      
      write(20000,*) av_div/(Nstep-Nstepskip), &
                    sqrt( (av_sq_div/(Nstep-Nstepskip))-(av_div/(Nstep-Nstepskip))**2)
      end program main

!*************************************************************************************!

!INITIALIZE
!*************************************************************************************! 
      subroutine INITIALIZE
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      Nstep  = 15500
      Nstepskip =5500
      Npolymer = 100
      Nbead_per_resi = 40
      Nresi = Npolymer*Nbead_per_resi
      Ndim           =3
      boxlen = 302.38
      Nbin = 500
      binsize = boxlen/(Nbin)
      rho     = Nresi/(boxlen**3)
      allocate(x(Nresi), y(Nresi), z(Nresi))
      allocate(com_lys(Npolymer, Ndim), com_glu(Npolymer,Ndim),com(Npolymer, Ndim))
      allocate(com_x(Npolymer,(Nstep-Nstepskip)))
      allocate(com_y(Npolymer,(Nstep-Nstepskip)))
      allocate(com_z(Npolymer,(Nstep-Nstepskip)))
      allocate(unnorm_summed_numbers(Nbin+1))
      unnorm_summed_numbers = 0
      av_droplet= 0
      av_div = 0
      av_sq_div=0
      allocate(num_neighbor(Npolymer),neighborlist(Npolymer,(Npolymer-1)))
      allocate(state(Npolymer))
      allocate(cluster_head(Npolymer),cluster_tail(Npolymer),cluster_size(Npolymer), cluster_checked(Npolymer))
      allocate (cluster_list(Npolymer))
      allocate(listed_id(Npolymer,Npolymer), id_droplet(Npolymer, Nstep))
      allocate(id_free_single(Npolymer, Nstep))
      allocate(com_droplet_x(Npolymer, Nstep))
      allocate(com_droplet_y(Npolymer, Nstep))
      allocate(com_droplet_z(Npolymer, Nstep))
      allocate(id_particle_in_cluster(Npolymer,Nstep))

      open (unit = 1 ,file &
       ="temp", &
                                                 action="read")
      !open(unit=900, file="cluster_dyn_0.3.dat", action = "write")
      !open(unit=1500, file="droplet_poly_1.0.dat", action = "write")     
      !open(unit=5000, file="msd_polymer_droplet.dat", action="write")
       open(unit=10000, file="anisotropy.dat",&
                                           action="write")
       open(unit=20000, file="anisotropy_averaged.dat",&
                                           action="write")
      end subroutine INITIALIZE

!*************************************************************************************!

!READ TRAJECTORY
!*************************************************************************************! 

      subroutine READ_FILE
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::  iresi
         read(1,*)
         read(1,"(i12)") step_num
         do iresi = 1, Nresi
            read(1,FMT1) x(iresi),y(iresi),z(iresi)
         end do
      !print*, x(Nresi),y (Nresi),z(Nresi)
      end subroutine READ_FILE


!*************************************************************************************!

!FINDING COM FOR POSITIVE AND NEGATIVE RESIDUES OF EACH CHAIN
!*************************************************************************************! 

      subroutine FINDCOM
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer::ii, ip,ll, iresi
      do ip =1 , Npolymer
         
         do ii = 1,3
            com_lys(ip,ii) =0
            com_glu(ip,ii) =0
            com(ip, ii)    =0
         end do
      end do
      do ip = 1, Npolymer
         ll = (ip-1)*Nbead_per_resi 
         do iresi = ll+1 , ll+10 
            com_lys(ip,1)= com_lys(ip,1)+x(iresi) 
            com_lys(ip,2)= com_lys(ip,2)+y(iresi)
            com_lys(ip,3)= com_lys(ip,3)+z(iresi)
         end do
         do iresi = ll +11, ll+20
            com_glu(ip,1) = com_glu(ip,1)+x(iresi)
            com_glu(ip,2) = com_glu(ip,2)+y(iresi)
            com_glu(ip,3) = com_glu(ip,3)+z(iresi)
         end do
      end do
      do ip = 1, Npolymer
         ll = (ip-1)*Nbead_per_resi
         do iresi = ll+1 , ll+Nbead_per_resi
            com(ip,1)= com(ip,1)+x(iresi)
            com(ip,2)= com(ip,2)+y(iresi)
            com(ip,3)= com(ip,3)+z(iresi)
         end do
      end do
      do ip = 1, Npolymer
         do  ii = 1,3
             com(ip,ii) =com(ip,ii)/Nbead_per_resi
             com_lys(ip,ii)=(2*com_lys(ip,ii)/Nbead_per_resi)
             com_glu(ip,ii) = (2*com_glu(ip,ii)/Nbead_per_resi)
         end do
      end do
      do ip =1, Npolymer
         com_x(ip,istep)=com(ip,1)
         com_y(ip, istep) = com(ip,2)
         com_z(ip, istep) = com(ip,3)
      end do
      end subroutine FINDCOM

!*************************************************************************************!

!BEAD to BEAD RDF
!*************************************************************************************! 

      subroutine RDF
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      real(kind=double) ::d(3),mod_d, pi=3.14, constant
      integer::iresi,jresi, whichbin
      integer:: ip,jp,ll,mm
      constant = (4/3)* pi
      do ip = 1, Npolymer-1
       do jp = ip+1, Npolymer
          if (ip .eq. jp) cycle
                d(1) = com(jp,1) -com(ip,1)
                d(2) = com(jp,2) -com(ip,2)
                d(3) = com(jp,3) -com(ip,3)
                mod_d = sqrt(d(1)**2 +d(2)**2 +d(3)**2)
                if (mod_d .gt. 3.7 .and. mod_d .lt. boxlen) then
                   !if(mod_d .lt. 3.8) then
                   !  write(500,*) iresi, jresi, mod_d
                   !end if
                   whichbin = int(mod_d/binsize)
                   unnorm_summed_numbers(whichbin+1) = unnorm_summed_numbers(whichbin+1)+2
                end if

       end do
      end do 
      end subroutine RDF
      
!*************************************************************************************!

!CLUSTERING TO SEE NEXT
!*************************************************************************************!

      subroutine GETCLUSTER
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      integer :: max_neighbor 
      integer :: ip, jp, pp, kk, mm, nbrs, nn,zz,qq,ll,rr
      integer :: ic,jc,ii
      real(kind=double) ::comij_x,comij_y, comij_z, modrij 
      real(kind=double):: fraction_droplet
      logical:: finding_not_completed
      integer :: temp, largest_cluster, droplet_polymer
      !Initializing at every time step
      clustersum = 0
      cluster_size = 0
      cluster_head  =0
      cluster_tail =0
      cluster_list =0
      cluster_checked =0 
      droplet_polymer=0
      do ip = 1, Npolymer
         do jp = 1, Npolymer-1
            neighborlist(ip, jp) = 0
         end do
      end do
      do ip = 1, Npolymer
         state(ip)        = 0
         num_neighbor(ip) = 0
      end do
      do ip =1, Npolymer
       do jp =1, Npolymer
          listed_id(ip,jp) =0
       end do
      end do
      
      !way to find the clusters
      !create the neighbor list 
      !go along the neighbor list to find all the neighbors of neighbors
      !which have not been counted to any cluster till now
      max_neighbor= Npolymer

      !search the neighbors
      do ip = 1, Npolymer -1
         do jp = ip+1, Npolymer
            comij_x = com(jp,1) -com(ip,1)
            comij_y = com(jp,2) -com(ip,2)
            comij_z = com(jp,3) -com(ip,3)
            modrij = sqrt (comij_x**2 + comij_y**2 + comij_z**2)
            if (modrij .le. 30) then 
               num_neighbor(ip) = num_neighbor(ip)+1
               num_neighbor(jp) = num_neighbor(jp)+1
               neighborlist(ip, num_neighbor(ip)) = jp
               neighborlist(jp, num_neighbor(jp)) = ip
            end if 
         end do
      end do
      
       
      
      !you have made your neighbor list 
      !now go on to find the clusters
      do ip = 1, Npolymer
         if (state(ip) .eq. 0) then
            clustersum =clustersum+1
            pp = ip
            state(pp) = clustersum
            cluster_head(clustersum) = pp
            cluster_tail(clustersum) = pp
            cluster_size(clustersum) = 1
            listed_id(clustersum, cluster_size(clustersum)) =pp
            id_particle_in_cluster(ip,istep) =clustersum
            !defined a logical variable in order to understand when to
            !end and when to start
            finding_not_completed  = .TRUE.
            do while (finding_not_completed)
               nbrs = num_neighbor(pp)
               kk = 1
               do while (kk .le. nbrs)
                  mm = neighborlist(pp, kk)
                  if (state(mm) .eq. 0) then 
                     cluster_list(cluster_tail(clustersum))= mm
                     cluster_tail(clustersum)              = mm
                     cluster_size(clustersum)              =cluster_size(clustersum)+1
                     state(mm) = clustersum
                     listed_id(clustersum, cluster_size(clustersum)) =mm
                     id_particle_in_cluster(mm,istep) = clustersum
                     !print*, listed_id(clustersum,cluster_size(clustersum))
                  end if
                  kk =kk+1
               end do
               cluster_checked(clustersum) = pp
               if (pp .eq.cluster_tail(clustersum)) then
                  finding_not_completed = .FALSE.
                  cluster_list(pp) = 0
               end if
               pp = cluster_list(pp)
            end do
         end if
      end do
      !print*, com(100,1)
      !do ic = 1, clustersum
      !   write(1100, *) cluster_size(ic)
      !end do
      do ic = 1, clustersum
         if(cluster_size(ic) .gt. 3) then
           droplet_polymer = droplet_polymer +cluster_size(ic)
         end if
      end do
      av_droplet= av_droplet+droplet_polymer
      fraction_droplet = droplet_polymer/Npolymer 
      largest_cluster= maxval(cluster_size) 
       do ic = 1, clustersum
         if( cluster_size(ic) .eq. largest_cluster) then
           zz = cluster_size(ic)
           !print*, zz
           do nn =1, zz
              id_droplet(listed_id(ic,nn),istep) =1
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)
              
           end do
         else 
          zz=cluster_size(ic)
          do nn =1, zz
              id_droplet(listed_id(ic,nn),istep) =0
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)
           end do
         end if
       end do

       do ic = 1, clustersum
         if( cluster_size(ic) .lt. 2) then
           zz = cluster_size(ic)
           !print*, zz
           do nn =1, zz
              id_free_single(listed_id(ic,nn),istep) =1
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)

           end do
         else 
           zz = cluster_size(ic)
           do nn =1, zz
              id_free_single(listed_id(ic,nn),istep) =0
              !print*,listed_id(ic,nn),id_droplet(listed_id(ic,nn),istep)

           end do
         end if
       end do





       !calculate the COM of droplet
       do ic = 1, clustersum

          do nn=1,cluster_size(ic)
             ll=listed_id(ic,nn)
             !qq=(ll-1)*Nbead_per_resi
             !do rr=1, Nbead_per_resi

              !  com_droplet_x(ic,istep)=com_droplet_x(ic,istep)+x(qq+rr)
              !  com_droplet_y(ic,istep)=com_droplet_y(ic,istep)+y(qq+rr)
              !  com_droplet_z(ic,istep)=com_droplet_z(ic,istep)+z(qq+rr)
               com_droplet_x(ic,istep)=com_droplet_x(ic,istep)+com_x(listed_id(ic,nn),istep)
               com_droplet_y(ic,istep)=com_droplet_y(ic,istep)+com_y(listed_id(ic,nn),istep)
               com_droplet_z(ic,istep)=com_droplet_z(ic,istep)+com_z(listed_id(ic,nn),istep)
             !end do
          end do
          !com_droplet_x(ic,istep)=com_droplet_x(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          !com_droplet_y(ic,istep)=com_droplet_y(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          !com_droplet_z(ic,istep)=com_droplet_z(ic,istep)/(Nbead_per_resi*cluster_size(ic))
          com_droplet_x(ic,istep)=com_droplet_x(ic,istep)/(cluster_size(ic))
          com_droplet_y(ic,istep)=com_droplet_y(ic,istep)/(cluster_size(ic))
          com_droplet_z(ic,istep)=com_droplet_z(ic,istep)/(cluster_size(ic))

       end do
 



       !Now find out the largest cluster

     


     !  do ic = 1, clustersum
     !    do jc = ic, clustersum 
     !       if(cluster_size(ic) .lt. cluster_size(jc)) then
     !         temp = cluster_size(jc)
     !         cluster_size(jc) = cluster_size(ic)
     !         cluster_size(ic) = temp
     !       end if
     !    end do
     ! end do
     ! largest_cluster = cluster_size(1)        
      !write(700,*) largest_cluster, clustersum
      !write(700,*) "stepdone=", istep 
      !write(800,*) step_num,real(Npolymer-droplet_polymer)/Npolymer,real(droplet_polymer)/Npolymer
      !write(900,*) step_num,largest_cluster

             
           
         

      end subroutine GETCLUSTER

!*************************************************************************************!

!ANISOTROPY OF THE DROPLET
!*************************************************************************************!
      subroutine GETANISOTROPY
      Use Precision_Control
      Use Param_Control
      Use Var_Control
      implicit none
      
      integer :: ip 
      real(kind=double):: d1,d2,d3,d_x_plus(Npolymer),d_y_plus(Npolymer), d_z_plus(Npolymer)
      real(kind=double):: d_x_minus(Npolymer),d_y_minus(Npolymer),d_z_minus(Npolymer)
      real(kind=double):: x_dim_max, y_dim_max, z_dim_max   
      real(kind=double):: x_dim_min, y_dim_min, z_dim_min   
      real(kind=double) :: tot_x, tot_y, tot_z     
      real(kind=double):: div1(3) 
      d_x_plus =0
      d_y_plus = 0
      d_z_plus =0
      d_x_minus=0
      d_y_minus=0
      d_z_minus=0
      div = 0


 
      do ip = 1, Npolymer
         if (id_droplet(ip,istep).eq. 1 ) then
             
             d1 =com_x(ip,istep)-com_droplet_x(id_particle_in_cluster(ip,istep), istep )
             d2=com_y(ip,istep)-com_droplet_y(id_particle_in_cluster(ip,istep), istep )
             d3=com_z(ip,istep)-com_droplet_z(id_particle_in_cluster(ip,istep), istep )
             if(d1 .ge. 0) then 
               d_x_plus(ip) = d1
             else if (d1 .lt. 0) then
               d_x_minus(ip) =abs(d1)
             end if


             if(d2 .ge. 0) then
               d_y_plus(ip) = d2
             else if (d2 .lt. 0) then
               d_y_minus(ip) =abs(d2)
             end if
             


             if(d3 .ge. 0) then
               d_z_plus(ip) = d3
             else if (d3 .lt. 0) then
               d_z_minus(ip) =abs(d3)
             end if
             
         end if 
      end do
      x_dim_max = maxval(d_x_plus)
      x_dim_min = -(maxval(d_x_minus))
      y_dim_max = maxval(d_y_plus)
      y_dim_min = -(maxval(d_y_minus))
      z_dim_max = maxval(d_z_plus)
      z_dim_min = -(maxval(d_z_minus))  

      tot_x= x_dim_max-x_dim_min
      tot_y = y_dim_max-y_dim_min
      tot_z = z_dim_max-z_dim_min
      
      div1(1) = tot_x/tot_y
      if (div1(1) .lt. 1.0) then
         div1(1) =1/div1(1)
      end if
   
      div1(2) = tot_y/tot_z
      if (div1(2) .lt. 1.0) then
         div1(2) =1/div1(2)
      end if

      div1(3) = tot_z/tot_x
      if (div1(3) .lt. 1.0) then
         div1(3) =1/div1(3)
      end if
      !div = maxval(div1)
       div = div1(1) +div1(2) +div1(3)
      !div= sqrt((tot_x-tot_y)**2 + (tot_x-tot_z)**2  + (tot_y-tot_z)**2)
      write(10000,*) istep, div
      
      end subroutine GETANISOTROPY
