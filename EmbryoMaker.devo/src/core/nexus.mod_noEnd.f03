!    EmbryoMaker software (General Node Model)
!    Computational model to simulate morphogenetic processes in living organs and tissues.
!    Copyright (C) 2014 Miquel Marin-Riera, Miguel Brun-Usan, Roland Zimm, Tommi Välikangas & Isaac Salazar-Ciudad

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.




module nexus

  use general
  use genetic
  use growth  !>>>> Miquel 17-6-13
  use death    !>>>> Miquel 18-6-13
  use ecm !>>>> Is 31-8-13 
  use single_node
  use mitosis
  use fitmo

contains

!******************************************************************

subroutine nexe

  integer i,j,k,ii,jj,kk,ik,ikk,ick,ret,othernode
  real*8 a,fitelli,refitGeneN,refitAngleN,req_i,req_j
  real*8 differ,pass
  real*8 kplast,kvol,dvol
  character*300 cx
  character*140 nofifit
  character*3 notwin 




  integer pid
  real*8 mz3,mz4,ai,bi
  character*1 q1
  character*2 q2
  character*3 q3
  character*4 q4
  character*5 q5
  character*6 q6
  real*8 :: maxg(15), mxa, mxb, mxxa, mxxb, ktime=0d0
  real*8,allocatable :: positions(:,:)
  character*140 :: thisfile, orifile, altrefile
  real*8 :: bbb,sco1,sco2
  real*8 :: scui,maxa,maxs(3),mins(3),u_d,u_e,minda
  integer :: nd_p,nd_a,nd_t,rtime_c,rtime_a,rtime_t,ndm,randl,auk,first=0
  real*8,allocatable :: mxg(:,:),ming(:,:)
  integer,allocatable :: dotpb(:)
  !   NEW AND TO TRy
  real*8 :: maxbasal(10),minbasal(10),maxapical(10),minapical(10),avgbasal(10),avgapical(10),nea,neb 
  character*20 :: whichend
!  integer,parameter :: icr=10
  integer,save::truenodeos=0
  real*8 :: icr,icrr=0.005 !RZ 30-05-17

  whichend="Real time 1"
  
!  print*,"now nexus"
 ! goto 171
  
! RZ 15-06-17: Let's save the true initial nodepars
if(truenodeos==0)then
truenodeos=1
nodeoini=0d0
  nodeoini(5,1)=node(3)%req ! epithelial
  nodeoini(5,2)=node(nd)%req ! mesenchyme
  nodeoini(6,1)=node(3)%da ! epithelial
  nodeoini(6,2)=node(nd)%da ! mesenchyme
  nodeoini(7,1)=node(3)%you ! epithelial
  nodeoini(7,2)=node(nd)%you ! mesenchyme
  nodeoini(8,1)=node(3)%adh ! epithelial
  nodeoini(8,2)=node(nd)%adh ! mesenchyme
  nodeoini(9,1)=node(3)%rep ! epithelial
  nodeoini(9,2)=node(nd)%rep ! mesenchyme
  nodeoini(10,1)=node(3)%req ! epithelial
  nodeoini(10,2)=node(nd)%req ! mesenchyme
  nodeoini(11,1)=node(3)%repcel ! epithelial
  nodeoini(11,2)=node(nd)%repcel ! mesenchyme
  nodeoini(12,1)=node(3)%tor ! epithelial
  nodeoini(12,2)=node(nd)%tor ! mesenchyme
  nodeoini(13,1)=node(3)%stor ! epithelial
  nodeoini(13,2)=node(nd)%stor ! mesenchyme
  nodeoini(14,1)=node(3)%reqs ! epithelial
  nodeoini(14,2)=node(nd)%reqs ! mesenchyme
  nodeoini(15,1)=node(3)%ke ! epithelial
  nodeoini(15,2)=node(nd)%ke ! mesenchyme
  nodeoini(16,1)=node(3)%mo ! epithelial
  nodeoini(16,2)=node(nd)%mo ! mesenchyme
  nodeoini(17,1)=node(3)%dmo ! epithelial
  nodeoini(17,2)=node(nd)%dmo ! mesenchyme
  nodeoini(28,1)=node(3)%kvol ! epithelial
  nodeoini(28,2)=node(nd)%kvol ! mesenchyme
  nodeoini(27,1)=node(3)%kplast ! epithelial
  nodeoini(27,2)=node(nd)%kplast ! mesenchyme
endif


! GENETIC REGULATION OF CELL BEHAVIOURS

  ! extracellular matrix secretion

  if (npag(nparam_per_node+4)>0) then ; call should_I_secrete   ; end if

  if(ffu(1)==0)then
    ! cell polarization  THIS ONE SHOULD BE THE FIRST IN HERE
    if (npag(nparam_per_node+8)>0) then ; call polarization ; end if
    !nparam_per_node+9 is to tell that growth is polarized 
    ! cell growth
    if (npag(nparam_per_node+1)>0) then ; call should_I_grow        ; end if  !It considers also polar growth in it, with nparam_per_node+9

    !cell division
    !if (npag(nparam_per_node+2)>0) then ; 
    call should_I_divide ; 
                                               !>>> Is 5-2-14

    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if

    ! change the size of the cell required for division
    if (npag(nparam_per_node+10)>0) then; call change_minsize_for_div ; end if

    !nparam_per_node+11 is to orient division according to the chemical polarization and not according to the physical (hertwig) one

    !nparam_per_node+12 the larger the most asymetric (in mass) is the plane of division

    ! change the maximal number of nodes per cell before the cell divides
    if (npag(nparam_per_node+15)>0) then; call change_maxsize_for_div ; end if  !>>> Is 23-3-14

  else
    if (ffu(25)==1) then ; if (npag(nparam_per_node+1)>0) then ; call should_I_grow ; end if  ; end if  !>>> Is 18-4-15 !It considers also polar growth in it, with nparam_per_node+9 ACHTUNG
    if (npag(nparam_per_node+8)>0) then ; call polarization_single ; end if
    if (npag(nparam_per_node+2)>0) then ; call should_I_divide_single ; end if
    !if (npag(nparam_per_node+8)>0) then ; call polarization_single ; end if
!print*,"req 2 nnn",node(nd-1)%req,nd-1
!print*,"req 1 nnn",node(nd)%req,nd
    ! cell apoptosis
    if (npag(nparam_per_node+3)>0.or.npag(nparam_per_node+14)>0) then ; call should_I_die  ; end if  !>>Miquel20-3-14
    if (npag(nparam_per_node+13)>0) then; call emt ; end if   ! >>> Is 4-10-14
   !call emt_single
  end if

  ! epithelial-mesenchymal transition
  if (npag(nparam_per_node+13)>0) then; call emt ; end if

  if (ffu(5)==1) then ; agex(:nd,1)=abs(node(:nd)%x) ;end if; !external source trick >>> Is 29-6-14
  if (ffu(5)==2) then ; agex(:nd,1)=node(:nd)%z**4    ;end if!((1.1d0+maxval(node(:nd)%z)+node(:nd)%z)**4) !external source trick miguel >>> Is 29-6-14

  if(ffu(5)==3)then !horizontal gradient (x-y plane) centered on the origin of coordinates
    do i=1,nd
      agex(i,6)=1/(1+(node(i)%x**2+node(i)%y**2)) !; print*,"i",i,"gex",agex(i,5)
    end do
  end if  
  
  if (ffu(15)==1)then; !especial conditions, gene 1 is always expressing in the borders !>>Miquel12-5-14
    do i=1,nd
      if(node(i)%hold==1)then
        agex(i,1)=1d0
        !if(node(i)%tipus<3)then
        !  agex(i,1)=1d0
        !else
        !  agex(i,4)=1d0
        !end if
      end if
    end do
  end if

  if (ffu(8)==1.and.nd>2) then   ! IS 23-4-13 this eliminates the nodes that get alone for too long
    ik=1
    do while(ik<=nd) !;print*,"node(",ik,")%talone=",node(ik)%talone
      if (node(ik)%talone>ttalone) then !;print*,ik,"entra mort",node(ik)%tipus
        if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
          ikk=node(ik)%icel !;print*,"entra mesenq"
          call apoptosis(ik)
          if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
            if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
          end if
        else
          if (node(node(ik)%altre)%talone>ttalone) then
            ikk=node(ik)%icel
            call apoptosis(ik)
            if (ikk>0) then   !notice that if ikk is negative it means that is not a celular node
              if(cels(ikk)%nunodes==0) call eliminate_cell(ikk)  !we eliminate the cell if it has no nodes left
            end if
          end if
        end if
        ik=ik+1 ! Is it? >>> Is 16-1-14
      else
        ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
      end if
    end do
!    call cellbreak ! to see if the cell is split in two 
  end if
  
  ! Check total number of nodes
 if(nd>20000)then; print*, "HAS BECOME TOO BIG"; 
  whichend="Over 20000 nd 3"
  goto 171
 endif 
  
  

! GENETIC REGULATION OF NODE PROPERTIES
  ! GENETIC REGULATION OF NODE PROPERTIES
  do i=1,nd         ! we update that parameter in each cell that expresses the gene
!print*, i, node(i)%tipus, node(i)%hold
                    ! WE ONLY UPDATES DE NODES IN WHICH THE GENE IS EXPRESSED, OTHERWISE WE LEAVE IT IS AS IT WAS 
    if(node(i)%hold==1) cycle  !we don't want the border cells to perform behaviours because that would alter and possibly break the border  !>>>>Miquel9-1-14
    if(node(i)%hold==2) cycle
    differ=1-node(i)%diffe !;print*,"differ",differ

    ! DIFFERENTIATION
    if (npag(25)>0) then
      ii=25;a=0.0; 
      do k=1,npag(ii) 
        kk=whonpag(ii,k) ; 
        if (gex(i,kk)>0.0d0) then ; 
          a=a+gex(i,kk)*gen(kk)%wa(ii)  
        endif
      end do
      ! NEW INDEPENDENT DIFFE !RZ 22-3-16
    if(ng.ge.11)  a=gen(11)%wa(25)*0.1 ! This makes that differ=0 in 10.000-100.000 iterations

      node(i)%diffe=node(i)%diffe+a*delta
      if (node(i)%diffe>1.0) node(i)%diffe=1.0 
      if (node(i)%diffe<0.0d0) node(i)%diffe=0.0

      nodeo(i)%diffe=node(i)%diffe !>>Miquel17-9-14
    end if
   ! print*, "t", whonpag(21,1)

    if(node(i)%tipus<2)then  !this for epithelial nodes  !>>Miquel12-5-14
      j=node(i)%altre

      if(ffu(11)==1) then        !plastic deformation   !>>Miquel5-2-14
        !kplast=node(i)%kplast
        !if(fmeanl(i)<epsilod)then;ki=1/node(i)%rep;else;ki=1/node(i)%you;endif
        !if(fmeanl(j)<epsilod)then;kj=1/node(j)%rep;else;kj=1/node(j)%you;endif

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aa=node(i)%reqp; bb=node(j)%reqp  ! previous values RZ 30-5-17

        node(i)%reqp=node(i)%reqp+(node(i)%kplast*fmeanl(i))*delta
        node(j)%reqp=node(j)%reqp+(node(j)%kplast*fmeanl(j))*delta !; print*,"fmeanl",fmeanl(i),fmeanl(j)

       ! RZ gradual 30-5-17
        if (node(i)%reqp>aa*(1+icr).and.aa.ne.0) node(i)%reqp=aa*(1+icr)
        if (node(j)%reqp>bb*(1+icr).and.bb.ne.0) node(j)%reqp=bb*(1+icr)
        if (node(i)%reqp<aa*(1-icr).and.aa.ne.0) node(i)%reqp=aa*(1-icr)
        if (node(j)%reqp<bb*(1-icr).and.bb.ne.0) node(j)%reqp=bb*(1-icr)
        if (aa==0.and.node(i)%reqp>icr) node(i)%reqp=icr
        if (bb==0.and.node(j)%reqp>icr) node(j)%reqp=icr
        if (aa==0.and.node(i)%reqp<icr*(-1)) node(i)%reqp=icr*(-1)
        if (bb==0.and.node(j)%reqp<icr*(-1)) node(j)%reqp=icr*(-1)
        ! RZ end changes

        nodeo(i)%reqp=node(i)%reqp  !>>Miquel17-9-14
        nodeo(j)%reqp=node(j)%reqp  !>>Miquel17-9-14

      else
        node(i)%reqp=0
        node(j)%reqp=0
      end if
    
      if (npag(21)>0) then  ! Contraction by genes 
        ii=21 ; a=0 ; b=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) !;print*,"ij",i,j,"gex",gex(i,kk),gex(j,kk),"kk",kk
          if (gex(i,kk)>0.0d0.or.gex(j,kk)>0.0d0) then 
            a=a+gex(i,kk)*gen(kk)%wa(ii)
            b=b+gex(j,kk)*gen(kk)%wa(ii)
          endif 
        enddo

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aaa=node(i)%reqc; bbb=node(j)%reqc  ! previous values RZ 30-5-17
        aa=0.25; bb=0.25!nodeoini(5,1); bb=nodeoini(5,1)  ! previous values RZ 30-5-17  
        ! GRAADUAL
        !if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa

        node(i)%reqc=nodeo(i)%reqc+a!*differ  !pfh no differ 26-11-15  
        node(j)%reqc=nodeo(j)%reqc+b!*differ  !pfh no differ 26-11-15 

        ! RZ gradual 30-5-17
        if (node(i)%reqc>aaa+aa*(icr).and.aa.ne.0) node(i)%reqc=aaa+aa*(icr)
        if (node(j)%reqc>bbb+bb*(icr).and.bb.ne.0) node(j)%reqc=bbb+bb*(icr)
        if (node(i)%reqc<aaa-aa*(icr).and.aa.ne.0) node(i)%reqc=aaa-aa*(icr)
        if (node(j)%reqc<bbb-bb*(icr).and.bb.ne.0) node(j)%reqc=bbb-bb*(icr)
        !if (aa==0.and.node(i)%reqc>delta*icr) node(i)%reqc=delta*icr
        !if (bb==0.and.node(j)%reqc>delta*icr) node(j)%reqc=delta*icr
        !if (aa==0.and.node(i)%reqc<delta*icr*-1) node(i)%reqc=delta*icr*-1
        !if (bb==0.and.node(j)%reqc<delta*icr*-1) node(j)%reqc=delta*icr*-1
        ! RZ end changes

      else
        node(i)%reqc=0 ;node(j)%reqc=0
      end if

      if(ffu(17)==1) then        !volume conservation   !>>Miquel6-5-14 ! >>> Is 24-5-14
        !kvol=node(i)%kvol
        dvol=0.5*(node(i)%reqcr+node(j)%reqcr-(node(i)%req+node(j)%req))

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aa=node(i)%reqv; bb=node(j)%reqv  ! previous values RZ 30-5-17

        node(i)%reqv=node(i)%reqv+node(i)%kvol*dvol*delta
        node(j)%reqv=node(j)%reqv+node(j)%kvol*dvol*delta


        ! RZ gradual 30-5-17
        if (node(i)%reqv>aa*(1+icr).and.aa.ne.0) node(i)%reqv=aa*(1+icr)
        if (node(j)%reqv>bb*(1+icr).and.bb.ne.0) node(j)%reqv=bb*(1+icr)
        if (node(i)%reqv<aa*(1-icr).and.aa.ne.0) node(i)%reqv=aa*(1-icr)
        if (node(j)%reqv<bb*(1-icr).and.bb.ne.0) node(j)%reqv=bb*(1-icr)
        if (aa==0.and.node(i)%reqv>icr) node(i)%reqv=icr
        if (bb==0.and.node(j)%reqv>icr) node(j)%reqv=icr
        if (aa==0.and.node(i)%reqv<icr*(-1)) node(i)%reqv=icr*(-1)
        if (bb==0.and.node(j)%reqv<icr*(-1)) node(j)%reqv=icr*(-1)
        ! RZ end changes

        nodeo(i)%reqv=node(i)%reqv  !>>Miquel17-9-14
        nodeo(j)%reqv=node(j)%reqv  !>>Miquel17-9-14

      else
        node(i)%reqv=0 ; node(j)%reqv=0
      end if

      if (ffu(18)==1) call diffusion_of_reqcr       !diffusion of reqcr ! >>> Is 25-5-14

      a=node(i)%da-node(i)%req !>>>Miquel9-9-15

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aa=node(i)%req; bb=node(j)%req  ! previous values RZ 30-5-17
      
      node(i)%req=node(i)%reqcr+node(i)%reqc+node(i)%reqp+node(i)%reqv  !now req is the sum of the req components: growth/apoptosis and contraction/deformation

!print*,"req 2 nnn222",node(nd-1)%req,nd-1
!print*,"req 1 nnn222",node(nd)%req,nd

      if(node(i)%req>df_reqmax) node(i)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%req<reqmin) node(i)%req=reqmin !the req can be deformed
      !node(i)%da=node(i)%req+a !>>>Miquel9-9-15
      nodeo(i)%da=node(i)%req+a !>>>Miquel9-9-15

      if(nodeo(i)%da.le.node(i)%req*1.0) nodeo(i)%da=node(i)%req*1.0  !pfh 9-6-15

      b=node(j)%da-node(j)%req !>>>Miquel9-9-15
      
      node(j)%req=node(j)%reqcr+node(j)%reqc+node(j)%reqp+node(j)%reqv  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(j)%req>df_reqmax) node(j)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(j)%req<reqmin) node(j)%req=reqmin !the req can be deformed
      !node(j)%da=node(j)%req+b !>>>Miquel9-9-15
      nodeo(j)%da=node(j)%req+b !>>>Miquel9-9-15

      if(node(i)%da.le.node(i)%req*1.0) node(i)%da=node(i)%req*1.0  !pfh 9-6-15 1.5
      if(node(i)%da.ge.node(i)%req*1.2) node(i)%da=node(i)%req*1.2  !pfh 9-6-15 2.5
      if(node(j)%da.le.node(j)%req*1.0) node(j)%da=node(j)%req*1.0  !pfh 9-6-15
      if(node(j)%da.ge.node(j)%req*1.2) node(j)%da=node(j)%req*1.2  !pfh 9-6-15


!pfh no differ 26-11-15  
      if (npag(6)>0) then
        ii=6;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
        if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
        if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%da) a=node(i)%da*0.1
      !if(b>0.1*node(j)%da) b=node(j)%da*0.1
      ! <<<18-5-17

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aaa=node(i)%da; bbb=node(j)%da  ! previous values RZ 30-5-17
        aa=nodeoini(6,1); bb=nodeoini(6,1)  ! previous values RZ 30-5-17


        !if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


        node(i)%da=nodeo(i)%da+a;node(j)%da=nodeo(j)%da+b

        ! RZ gradual 30-5-17
        if (node(i)%da>aaa+aa*(icr).and.aa.ne.0) node(i)%da=aaa+aa*(icr)
        if (node(j)%da>bbb+bb*(icr).and.bb.ne.0) node(j)%da=bbb+bb*(icr)
        if (node(i)%da<aaa-aa*(icr).and.aa.ne.0) node(i)%da=aaa-aa*(icr)
        if (node(j)%da<bbb-bb*(icr).and.bb.ne.0) node(j)%da=bbb-bb*(icr)
        !if (aa==0.and.node(i)%da>delta*icr) node(i)%da=delta*icr
        !if (bb==0.and.node(j)%da>delta*icr) node(j)%da=delta*icr
        !if (aa==0.and.node(i)%da<delta*icr*-1) node(i)%da=delta*icr*-1
        !if (bb==0.and.node(j)%da<delta*icr*-1) node(j)%da=delta*icr*-1
        ! RZ end changes


        if (node(i)%da<0) node(i)%da=0.0;if (node(j)%da<0) node(j)%da=0.0

      if(node(i)%da.le.node(i)%req*1.0) node(i)%da=node(i)%req*1.0  !pfh 9-6-15
      if(node(i)%da.ge.node(i)%req*1.2) node(i)%da=node(i)%req*1.2  !pfh 9-6-15
      if(node(j)%da.le.node(j)%req*1.0) node(j)%da=node(j)%req*1.0  !pfh 9-6-15
      if(node(j)%da.ge.node(j)%req*1.2) node(j)%da=node(j)%req*1.2  !pfh 9-6-15

      else

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%da) a=node(i)%da*0.1
      !if(b>0.1*node(j)%da) b=node(j)%da*0.1
      ! <<<18-5-17

        ! Rz 18-5-17 makes node properties to change gradually?????
         icr=icrr ! RZ 30-5-17
        aaa=node(i)%da; bbb=node(j)%da  ! previous values RZ 30-5-17
        aa=nodeoini(6,1); bb=nodeoini(6,1)  ! previous values RZ 30-5-17


        !if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


        node(i)%da=nodeo(i)%da ; node(j)%da=nodeo(j)%da !pfh 9-6-15

        ! RZ gradual 30-5-17
        if (node(i)%da>aaa+aa*(icr).and.aa.ne.0) node(i)%da=aaa+aa*(icr)
        if (node(j)%da>bbb+bb*(icr).and.bb.ne.0) node(j)%da=bbb+bb*(icr)
        if (node(i)%da<aaa-aa*(icr).and.aa.ne.0) node(i)%da=aaa-aa*(icr)
        if (node(j)%da<bbb-bb*(icr).and.bb.ne.0) node(j)%da=bbb-bb*(icr)
        !if (aa==0.and.node(i)%da>delta*icr) node(i)%da=delta*icr
        !if (bb==0.and.node(j)%da>delta*icr) node(j)%da=delta*icr
        !if (aa==0.and.node(i)%da<delta*icr*-1) node(i)%da=delta*icr*-1
        !if (bb==0.and.node(j)%da<delta*icr*-1) node(j)%da=delta*icr*-1
        ! RZ end changes

      end if

      if(node(i)%da.le.node(i)%req*1.0) node(i)%da=node(i)%req*1.0  !pfh 9-6-15 1.5
      if(node(i)%da.ge.node(i)%req*1.2) node(i)%da=node(i)%req*1.2  !pfh 9-6-15 2.5
      if(node(j)%da.le.node(j)%req*1.0) node(j)%da=node(j)%req*1.0  !pfh 9-6-15
      if(node(j)%da.ge.node(j)%req*1.2) node(j)%da=node(j)%req*1.2  !pfh 9-6-15

      if (npag(7)>0) then;ii=7;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  

      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif

      enddo

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%you; bbb=node(j)%you  ! previous values RZ 30-5-17
      aa=nodeoini(7,1); bb=nodeoini(7,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%you) a=node(i)%you*0.1
      !if(b>0.1*node(j)%you) b=node(j)%you*0.1
      ! <<<18-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%you=nodeo(i)%you+a;node(j)%you=nodeo(j)%you+b

       ! RZ gradual 30-5-17
        if (node(i)%you>aaa+aa*(icr).and.aa.ne.0) node(i)%you=aaa+aa*(icr)
        if (node(j)%you>bbb+bb*(icr).and.bb.ne.0) node(j)%you=bbb+bb*(icr)
        if (node(i)%you<aaa-aa*(icr).and.aa.ne.0) node(i)%you=aaa-aa*(icr)
        if (node(j)%you<bbb-bb*(icr).and.bb.ne.0) node(j)%you=bbb-bb*(icr)
        !if (aa==0.and.node(i)%you>delta*icr) node(i)%you=delta*icr
        !if (bb==0.and.node(j)%you>delta*icr) node(j)%you=delta*icr
        !if (aa==0.and.node(i)%you<delta*icr*-1) node(i)%you=delta*icr*-1
       ! if (bb==0.and.node(j)%you<delta*icr*-1) node(j)%you=delta*icr*-1
        ! RZ end changes

      if (node(i)%you<0) node(i)%you=0.0;if (node(j)%you<0) node(j)%you=0.0;end if

      if (npag(8)>0) then;ii=8;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;
if(3==30)then
         icr=icrr
        if (nodeo(i)%adh+a>nodeo(i)%adh*(1+delta*icr)) a=node(i)%adh*(delta*icr)
        if (nodeo(j)%adh+b>nodeo(j)%adh*(1+delta*icr)) b=node(j)%adh*(delta*icr)
        if (nodeo(i)%adh+a<nodeo(i)%adh*(1-delta*icr)) a=node(i)%adh*(delta*icr)*(-1)
        if (nodeo(j)%adh+b<nodeo(j)%adh*(1-delta*icr)) b=node(j)%adh*(delta*icr)*(-1)
        if (nodeo(i)%adh==0.and.a>delta*icr) a=delta*icr
        if (nodeo(j)%adh==0.and.b>delta*icr) b=delta*icr
        if (nodeo(i)%adh==0.and.a<-1*delta*icr) a=-1*delta*icr
        if (nodeo(j)%adh==0.and.b<-1*delta*icr) b=-1*delta*icr
endif
      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%adh) a=node(i)%adh*0.1
      !if(b>0.1*node(j)%adh) b=node(j)%adh*0.1
      ! <<<18-5-17

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%adh; bbb=node(j)%adh  ! previous values RZ 30-5-17
       aa=nodeoini(8,1); bb=nodeoini(8,1)  ! previous values RZ 30-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%adh=nodeo(i)%adh+a;node(j)%adh=nodeo(j)%adh+b

       ! RZ gradual 30-5-17
        if (node(i)%adh>aaa+aa*(icr).and.aa.ne.0) node(i)%adh=aaa+aa*(icr)
        if (node(j)%adh>bbb+bb*(icr).and.bb.ne.0) node(j)%adh=bbb+bb*(icr)
        if (node(i)%adh<aaa-aa*(icr).and.aa.ne.0) node(i)%adh=aaa-aa*(icr)
        if (node(j)%adh<bbb-bb*(icr).and.bb.ne.0) node(j)%adh=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%adh>delta*icr) node(i)%adh=delta*icr
       ! if (bb==0.and.node(j)%adh>delta*icr) node(j)%adh=delta*icr
       ! if (aa==0.and.node(i)%adh<delta*icr*-1) node(i)%adh=delta*icr*-1
       ! if (bb==0.and.node(j)%adh<delta*icr*-1) node(j)%adh=delta*icr*-1
        ! RZ end changes

      if (node(i)%adh<0) node(i)%adh=0.0;if (node(j)%adh<0) node(j)%adh=0.0;end if

      if (npag(9)>0) then;ii=9;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%rep; bbb=node(j)%rep  ! previous values RZ 30-5-17
      aa=nodeoini(9,1); bb=nodeoini(9,1)  ! previous values RZ 30-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%rep) a=node(i)%rep*0.1
      !if(b>0.1*node(j)%rep) b=node(j)%rep*0.1
      ! <<<18-5-17

      node(i)%rep=nodeo(i)%rep+a;node(j)%rep=nodeo(j)%rep+b

       ! RZ gradual 30-5-17
        if (node(i)%rep>aaa+aa*(icr).and.aa.ne.0) node(i)%rep=aaa+aa*(icr)
        if (node(j)%rep>bbb+bb*(icr).and.bb.ne.0) node(j)%rep=bbb+bb*(icr)
        if (node(i)%rep<aaa-aa*(icr).and.aa.ne.0) node(i)%rep=aaa-aa*(icr)
        if (node(j)%rep<bbb-bb*(icr).and.bb.ne.0) node(j)%rep=bbb-bb*(icr)
        !if (aa==0.and.node(i)%rep>delta*icr) node(i)%rep=delta*icr
        !if (bb==0.and.node(j)%rep>delta*icr) node(j)%rep=delta*icr
        !if (aa==0.and.node(i)%rep<delta*icr*-1) node(i)%rep=delta*icr*-1
        !if (bb==0.and.node(j)%rep<delta*icr*-1) node(j)%rep=delta*icr*-1
        ! RZ end changes

      if (node(i)%rep<0) node(i)%rep=0.0;if (node(j)%rep<0) node(j)%rep=0.0;end if

      if (npag(10)>0) then;ii=10;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%repcel; bbb=node(j)%repcel  ! previous values RZ 30-5-17
      aa=nodeoini(10,1); bb=nodeoini(10,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%repcel) a=node(i)%repcel*0.1
      !if(b>0.1*node(j)%repcel) b=node(j)%repcel*0.1
      ! <<<18-5-17
       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%repcel=nodeo(i)%repcel+a;node(j)%repcel=nodeo(j)%repcel+b

       ! RZ gradual 30-5-17
        if (node(i)%repcel>aaa+aa*(icr).and.aa.ne.0) node(i)%repcel=aaa+aa*(icr)
        if (node(j)%repcel>bbb+bb*(icr).and.bb.ne.0) node(j)%repcel=bbb+bb*(icr)
        if (node(i)%repcel<aaa-aa*(icr).and.aa.ne.0) node(i)%repcel=aaa-aa*(icr)
        if (node(j)%repcel<bbb-bb*(icr).and.bb.ne.0) node(j)%repcel=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%repcel>delta*icr) node(i)%repcel=delta*icr
       ! if (bb==0.and.node(j)%repcel>delta*icr) node(j)%repcel=delta*icr
       ! if (aa==0.and.node(i)%repcel<delta*icr*-1) node(i)%repcel=delta*icr*-1
       ! if (bb==0.and.node(j)%repcel<delta*icr*-1) node(j)%repcel=delta*icr*-1
        ! RZ end changes

      if (node(i)%repcel<0) node(i)%repcel=0.0;if (node(j)%repcel<0) node(j)%repcel=0.0;end if

      if (npag(11)>0) then;ii=11;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%tor; bbb=node(j)%tor  ! previous values RZ 30-5-17
      aa=nodeoini(11,1); bb=nodeoini(11,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%tor) a=node(i)%tor*0.1
      !if(b>0.1*node(j)%tor) b=node(j)%tor*0.1
      ! <<<18-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%tor=nodeo(i)%tor+a;node(j)%tor=nodeo(j)%tor+b

       ! RZ gradual 30-5-17
        if (node(i)%tor>aaa+aa*(icr).and.aa.ne.0) node(i)%tor=aaa+aa*(icr)
        if (node(j)%tor>bbb+bb*(icr).and.bb.ne.0) node(j)%tor=bbb+bb*(icr)
        if (node(i)%tor<aaa-aa*(icr).and.aa.ne.0) node(i)%tor=aaa-aa*(icr)
        if (node(j)%tor<bbb-bb*(icr).and.bb.ne.0) node(j)%tor=bbb-bb*(icr)
        !if (aa==0.and.node(i)%tor>delta*icr) node(i)%tor=delta*icr
        !if (bb==0.and.node(j)%tor>delta*icr) node(j)%tor=delta*icr
        !if (aa==0.and.node(i)%tor<delta*icr*-1) node(i)%tor=delta*icr*-1
        !if (bb==0.and.node(j)%tor<delta*icr*-1) node(j)%tor=delta*icr*-1
        ! RZ end changes

      if (node(i)%tor<0) node(i)%tor=0.0;if (node(j)%tor<0) node(j)%tor=0.0;end if

      if (npag(12)>0) then;ii=12;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%stor; bbb=node(j)%stor  ! previous values RZ 30-5-17
      aa=nodeoini(12,1); bb=nodeoini(12,1)  ! previous values RZ 30-5-17
      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%stor) a=node(i)%stor*0.1
      !if(b>0.1*node(j)%stor) b=node(j)%stor*0.1
      ! <<<18-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%stor=nodeo(i)%stor+a;node(j)%stor=nodeo(j)%stor+b

       ! RZ gradual 30-5-17
        if (node(i)%stor>aaa+aa*(icr).and.aa.ne.0) node(i)%stor=aaa+aa*(icr)
        if (node(j)%stor>bbb+bb*(icr).and.bb.ne.0) node(j)%stor=bbb+bb*(icr)
        if (node(i)%stor<aaa-aa*(icr).and.aa.ne.0) node(i)%stor=aaa-aa*(icr)
        if (node(j)%stor<bbb-bb*(icr).and.bb.ne.0) node(j)%stor=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%stor>delta*icr) node(i)%stor=delta*icr
       ! if (bb==0.and.node(j)%stor>delta*icr) node(j)%stor=delta*icr
       ! if (aa==0.and.node(i)%stor<delta*icr*-1) node(i)%stor=delta*icr*-1
       ! if (bb==0.and.node(j)%stor<delta*icr*-1) node(j)%stor=delta*icr*-1
        ! RZ end changes

      if (node(i)%stor<0) node(i)%stor=0.0;if (node(j)%stor<0) node(j)%stor=0.0;end if

      if (npag(13)>0) then;ii=13;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%reqs; bbb=node(j)%reqs  ! previous values RZ 30-5-17
      aa=nodeoini(13,1); bb=nodeoini(13,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%reqs) a=node(i)%reqs*0.1
      !if(b>0.1*node(j)%reqs) b=node(j)%reqs*0.1
      ! <<<18-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%reqs=nodeo(i)%reqs+a;node(j)%reqs=nodeo(j)%reqs+b

       ! RZ gradual 30-5-17
        if (node(i)%reqs>aaa+aa*(icr).and.aa.ne.0) node(i)%reqs=aaa+aa*(icr)
        if (node(j)%reqs>bbb+bb*(icr).and.bb.ne.0) node(j)%reqs=bbb+bb*(icr)
        if (node(i)%reqs<aaa-aa*(icr).and.aa.ne.0) node(i)%reqs=aaa-aa*(icr)
        if (node(j)%reqs<bbb-bb*(icr).and.bb.ne.0) node(j)%reqs=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%reqs>delta*icr) node(i)%reqs=delta*icr
       ! if (bb==0.and.node(j)%reqs>delta*icr) node(j)%reqs=delta*icr
       ! if (aa==0.and.node(i)%reqs<delta*icr*-1) node(i)%reqs=delta*icr*-1
       ! if (bb==0.and.node(j)%reqs<delta*icr*-1) node(j)%reqs=delta*icr*-1
        ! RZ end changes

      if (node(i)%reqs<0) node(i)%reqs=0.0;if (node(j)%reqs<0) node(j)%reqs=0.0;end if

      if (npag(14)>0) then;ii=14;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%ke; bbb=node(j)%ke  ! previous values RZ 30-5-17
      aa=nodeoini(14,1); bb=nodeoini(14,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%ke) a=node(i)%ke*0.1
      !if(b>0.1*node(j)%ke) b=node(j)%ke*0.1
      ! <<<18-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%ke=nodeo(i)%ke+a;node(j)%ke=nodeo(j)%ke+b

       ! RZ gradual 30-5-17
        if (node(i)%ke>aaa+aa*(icr).and.aa.ne.0) node(i)%ke=aaa+aa*(icr)
        if (node(j)%ke>bbb+bb*(icr).and.bb.ne.0) node(j)%ke=bbb+bb*(icr)
        if (node(i)%ke<aaa-aa*(icr).and.aa.ne.0) node(i)%ke=aaa-aa*(icr)
        if (node(j)%ke<bbb-bb*(icr).and.bb.ne.0) node(j)%ke=bbb-bb*(icr)
        !if (aa==0.and.node(i)%ke>delta*icr) node(i)%ke=delta*icr
        !if (bb==0.and.node(j)%ke>delta*icr) node(j)%ke=delta*icr
        !if (aa==0.and.node(i)%ke<delta*icr*-1) node(i)%ke=delta*icr*-1
        !if (bb==0.and.node(j)%ke<delta*icr*-1) node(j)%ke=delta*icr*-1
        ! RZ end changes

      if (node(i)%ke<0) node(i)%ke=0.0;if (node(j)%ke<0) node(j)%ke=0.0;end if

      if (npag(15)>0) then;ii=15;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in force units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%mo; bbb=node(j)%mo  ! previous values RZ 30-5-17
      aa=nodeoini(15,1); bb=nodeoini(15,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%mo) a=node(i)%mo*0.1
      !if(b>0.1*node(j)%mo) b=node(j)%mo*0.1
      ! <<<18-5-17
       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%mo=nodeo(i)%mo+a;node(j)%mo=nodeo(j)%mo+b

      ! RZ gradual 30-5-17
        if (node(i)%mo>aaa+aa*(icr).and.aa.ne.0) node(i)%mo=aaa+aa*(icr)
        if (node(j)%mo>bbb+bb*(icr).and.bb.ne.0) node(j)%mo=bbb+bb*(icr)
        if (node(i)%mo<aaa-aa*(icr).and.aa.ne.0) node(i)%mo=aaa-aa*(icr)
        if (node(j)%mo<bbb-bb*(icr).and.bb.ne.0) node(j)%mo=bbb-bb*(icr)
        !if (aa==0.and.node(i)%mo>delta*icr) node(i)%mo=delta*icr
        !if (bb==0.and.node(j)%mo>delta*icr) node(j)%mo=delta*icr
        !if (aa==0.and.node(i)%mo<delta*icr*-1) node(i)%mo=delta*icr*-1
        !if (bb==0.and.node(j)%mo<delta*icr*-1) node(j)%mo=delta*icr*-1
        ! RZ end changes

      if (node(i)%mo<0) node(i)%mo=0.0;if (node(j)%mo<0) node(j)%mo=0.0;
      else; if(node(i)%mo==10) node(i)%mo=1000 !!!!!!
      ! RZ NEW, JUST TO ENSURE
      end if

      if (npag(16)>0) then;ii=16;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif
      enddo

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%dmo; bbb=node(j)%dmo  ! previous values RZ 30-5-17
      aa=nodeoini(16,1); bb=nodeoini(16,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%dmo) a=node(i)%dmo*0.1
      !if(b>0.1*node(j)%dmo) b=node(j)%dmo*0.1
      ! <<<18-5-17
       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
       ! if(b>(delta*icr)*aa) b=(delta*icr)*aa
       ! if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%dmo=nodeo(i)%dmo+a;node(j)%dmo=nodeo(j)%dmo+b

    ! RZ gradual 30-5-17
        if (node(i)%dmo>aaa+aa*(icr).and.aa.ne.0) node(i)%dmo=aaa+aa*(icr)
        if (node(j)%dmo>bbb+bb*(icr).and.bb.ne.0) node(j)%dmo=bbb+bb*(icr)
        if (node(i)%dmo<aaa-aa*(icr).and.aa.ne.0) node(i)%dmo=aaa-aa*(icr)
        if (node(j)%dmo<bbb-bb*(icr).and.bb.ne.0) node(j)%dmo=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%dmo>delta*icr) node(i)%dmo=delta*icr
       ! if (bb==0.and.node(j)%dmo>delta*icr) node(j)%dmo=delta*icr
       ! if (aa==0.and.node(i)%dmo<delta*icr*-1) node(i)%dmo=delta*icr*-1
       ! if (bb==0.and.node(j)%dmo<delta*icr*-1) node(j)%dmo=delta*icr*-1
        ! RZ end changes

      if (node(i)%dmo<0) node(i)%dmo=0.0;if (node(j)%dmo<0) node(j)%dmo=0.0;end if

      if (npag(27)>0) then; ii=27;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units  
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif;enddo;

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%kplast; bbb=node(j)%kplast  ! previous values RZ 30-5-17
      aa=nodeoini(27,1); bb=nodeoini(27,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%kplast) a=node(i)%kplast*0.1
      !if(b>0.1*node(j)%kplast) b=node(j)%kplast*0.1
      ! <<<18-5-17

        !if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%kplast=nodeo(i)%kplast+a;node(j)%kplast=nodeo(j)%kplast+b

    ! RZ gradual 30-5-17
        if (node(i)%kplast>aaa+aa*(icr).and.aa.ne.0) node(i)%kplast=aaa+aa*(icr)
        if (node(j)%kplast>bbb+bb*(icr).and.bb.ne.0) node(j)%kplast=bbb+bb*(icr)
        if (node(i)%kplast<aaa-aa*(icr).and.aa.ne.0) node(i)%kplast=aaa-aa*(icr)
        if (node(j)%kplast<bbb-bb*(icr).and.bb.ne.0) node(j)%kplast=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%kplast>delta*icr) node(i)%kplast=delta*icr
       ! if (bb==0.and.node(j)%kplast>delta*icr) node(j)%kplast=delta*icr
       ! if (aa==0.and.node(i)%kplast<delta*icr*-1) node(i)%kplast=delta*icr*-1
       ! if (bb==0.and.node(j)%kplast<delta*icr*-1) node(j)%kplast=delta*icr*-1
        ! RZ end changes

      if (node(i)%kplast<0) node(i)%kplast=0.0;if (node(j)%kplast<0) node(j)%kplast=0.0;end if

      if (npag(28)>0) then;ii=28;a=0 ; b=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; !wa in req-space units 
      if (gex(i,kk)>0.0d0) then ; a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;
      if (gex(j,kk)>0.0d0) then ; b=b+gex(j,kk)*gen(kk)%wa(ii) ;endif
      enddo
 
      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%kvol; bbb=node(j)%kvol  ! previous values RZ 30-5-17
      aa=nodeoini(28,1); bb=nodeoini(28,1)  ! previous values RZ 30-5-17

      ! >>>18-5-17 gradual
      !if(a>0.1*node(i)%kvol) a=node(i)%kvol*0.1
      !if(b>0.1*node(j)%kvol) b=node(j)%kvol*0.1
      ! <<<18-5-17

        !if(a>(delta*icr)*aa) a=(delta*icr)*aa
        !if(a<(-delta*icr)*aa) a=(-delta*icr)*aa
        !if(b>(delta*icr)*aa) b=(delta*icr)*aa
        !if(b<(-delta*icr)*aa) b=(-delta*icr)*aa


      node(i)%kvol=nodeo(i)%kvol+a;node(j)%kvol=nodeo(j)%kvol+b

    ! RZ gradual 30-5-17
       if (node(i)%kvol>aaa+aa*(icr).and.aa.ne.0) node(i)%kvol=aaa+aa*(icr)
       if (node(j)%kvol>bbb+bb*(icr).and.bb.ne.0) node(j)%kvol=bbb+bb*(icr)
       if (node(i)%kvol<aaa-aa*(icr).and.aa.ne.0) node(i)%kvol=aaa-aa*(icr)
       if (node(j)%kvol<bbb-bb*(icr).and.bb.ne.0) node(j)%kvol=bbb-bb*(icr)
       ! if (aa==0.and.node(i)%kvol>delta*icr) node(i)%kvol=delta*icr
       ! if (bb==0.and.node(j)%kvol>delta*icr) node(j)%kvol=delta*icr
       ! if (aa==0.and.node(i)%kvol<delta*icr*-1) node(i)%kvol=delta*icr*-1
       ! if (bb==0.and.node(j)%kvol<delta*icr*-1) node(j)%kvol=delta*icr*-1
        ! RZ end changes


      if (node(i)%kvol<0) node(i)%kvol=0.0;if (node(j)%kvol<0) node(j)%kvol=0.0; end if

      ! NEW NEW NEW: We ensure the da is equal in both faces
     !  if(node(i)%req<
      reqmin=0.1 !w1
      screen_radius=0.85 !w!
      !if(node(i)%stor>10*node(i)%tor) node(i)%stor=10*node(i)%tor
      if(node(i)%tor>10*node(i)%stor) node(i)%tor=10*node(i)%stor
      if(node(i)%stor>500) node(i)%stor=500
      if(node(i)%tor>500) node(i)%tor=500
      if(node(i)%stor<1) node(i)%stor=1
      if(node(i)%tor<1) node(i)%tor=1
      !if(node(j)%stor>10*node(j)%tor) node(j)%stor=10*node(j)%tor
      if(node(j)%tor>10*node(j)%stor) node(j)%tor=10*node(j)%stor
      if(node(j)%stor>500) node(j)%stor=500
      if(node(j)%tor>500) node(j)%tor=500
      if(node(j)%stor<1) node(j)%stor=1
      if(node(j)%tor<1) node(j)%tor=1
     ! if(node(j)%stor>200) node(j)%stor=200
      !if(node(j)%tor>20) node(j)%tor=20
      !node(i)%stor=10; node(i)%tor=100
      !node(j)%stor=10; node(j)%tor=100
      !node(i)%stor=90d0; node(i)%tor=10d0 !w!
      !if(node(i)%da<node(j)%da) node(i)%da=node(j)%da
     ! if(node(j)%da<node(i)%da) node(j)%da=node(i)%da
    ! if(node(i)%stor<0.5*node(i)%tor) node(i)%stor=0.5*node(i)%tor
     ! if(node(i)%stor>2*node(i)%tor) node(i)%stor=2*node(i)%tor
      !print*, i, node(i)%tipus, node(i)%tor, node(i)%stor, node(i)%req, node(i)%da
      ! scp rolazimm@128.214.53.22:/home/rolazimm/Desktop/ensemble/originals/ori2_selected/en*ori2.gz ensi09ori2
      !node(i)%adh=10d0

!print*, i, "111", npag(21)

!print*,"req 2 nnn333",node(nd-1)%req,nd-1
!print*,"req 1 nnn333",node(nd)%req,nd

    else if (node(i)%tipus>2)then !this for mesenchyme and ECM ! NEW NEW NEW

!cycle ! OUTCOMMENT 

!print*, i, "222", npag(21)

      if (npag(21)>0) then 
        ii=21 ; a=0
        do k=1,npag(ii) ; 
          kk=whonpag(ii,k) 
          if (gex(i,kk)>0.0d0) then ; 
            a=a+gex(i,kk)*gen(kk)%wa(ii)

          endif 
        enddo

      ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%reqc
      aa=nodeoini(5,2)  ! previous values RZ 30-5-17

       ! if(a>(delta*icr)*aa) a=(delta*icr)*aa
       ! if(a<(-delta*icr)*aa) a=(-delta*icr)*aa

        !node(i)%reqc=nodeo(i)%reqc+a*differ*delta !wa in req-space units !>>>Miquel 9-9-15
        node(i)%reqc=nodeo(i)%reqc+a !>>>Miquel 9-9-15 !pfh no differ 26-11-15  

    ! RZ gradual 30-5-17
        if (node(i)%reqc>aaa+aa*(icr).and.aa.ne.0) node(i)%reqc=aaa+aa*(icr)
        if (node(i)%reqc<aaa-aa*(icr).and.aa.ne.0) node(i)%reqc=aaa-aa*(icr)
       ! if (aa==0.and.node(i)%reqc>delta*icr) node(i)%reqc=delta*icr
       ! if (aa==0.and.node(i)%reqc<delta*icr*-1) node(i)%reqc=delta*icr*-1
        ! RZ end changes

      else
        node(i)%reqc=0
      end if

      !a=node(i)%da-node(i)%req !>>>Miquel9-9-15
      a=nodeo(i)%da-node(i)%req !>>>Miquel9-9-15
      node(i)%req=node(i)%reqcr+node(i)%reqc  !now req is the sum of the req components: growth/apoptosis and contraction/deformation
      if(node(i)%req>df_reqmax) node(i)%req=df_reqmax !put an upper an lower boundary on how much  !>>Miquel28-7-14
      if(node(i)%req<reqmin) node(i)%req=reqmin !the req can be deformed !NEW *5 deleted
      !node(i)%da=node(i)%req+a !>>>Miquel9-9-15

      nodeo(i)%da=node(i)%req+a !>>>Miquel9-9-15

      if(nodeo(i)%da.le.node(i)%req*1.0) nodeo(i)%da=node(i)%req*1.0  !pfh 9-6-15
      if(nodeo(i)%da.ge.node(i)%req*1.2) nodeo(i)%da=node(i)%req*1.2  !pfh 9-6-15

!pfh no differ 26-11-15  
      if (npag(6)>0) then;ii=6;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
        a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%da
       aa=nodeoini(6,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%da=nodeo(i)%da+a

    ! RZ gradual 30-5-17
        if (node(i)%da>aaa+aa*(icr).and.aa.ne.0) node(i)%da=aaa+aa*(icr)
        if (node(i)%da<aaa-aa*(icr).and.aa.ne.0) node(i)%da=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%da>delta*icr) node(i)%da=delta*icr
    !    if (aa==0.and.node(i)%da<delta*icr*-1) node(i)%da=delta*icr*-1
    ! RZ end changes

      if (node(i)%da<0) node(i)%da=0.0
      if(node(i)%da.le.node(i)%req*1.0) node(i)%da=node(i)%req*1.0  !pfh 9-6-15
      if(node(i)%da.ge.node(i)%req*1.2) node(i)%da=node(i)%req*1.2  !pfh 9-6-15
      else
        node(i)%da=nodeo(i)%da !pfh 9-6-15
      end if
      if (npag(7)>0) then;ii=7;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%you
      aa=nodeoini(7,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%you=nodeo(i)%you+a

    ! RZ gradual 30-5-17
        if (node(i)%you>aaa+aa*(icr).and.aa.ne.0) node(i)%you=aaa+aa*(icr)
        if (node(i)%you<aaa-aa*(icr).and.aa.ne.0) node(i)%you=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%you>delta*icr) node(i)%you=delta*icr
    !    if (aa==0.and.node(i)%you<delta*icr*-1) node(i)%you=delta*icr*-1
    ! RZ end changes

      if (node(i)%you<0) node(i)%you=0.0;end if
      if (npag(8)>0) then;ii=8;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%adh
       aa=nodeoini(8,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
    !    if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%adh=nodeo(i)%adh+a

    ! RZ gradual 30-5-17
        if (node(i)%adh>aaa+aa*(icr).and.aa.ne.0) node(i)%adh=aaa+aa*(icr)
        if (node(i)%adh<aaa-aa*(icr).and.aa.ne.0) node(i)%adh=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%adh>delta*icr) node(i)%adh=delta*icr
    !    if (aa==0.and.node(i)%adh<delta*icr*-1) node(i)%adh=delta*icr*-1
    ! RZ end changes

      if (node(i)%adh<0) node(i)%adh=0.0;end if
      if (npag(9)>0) then;ii=9;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%rep
      aa=nodeoini(9,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
    !    if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%rep=nodeo(i)%rep+a

    ! RZ gradual 30-5-17
        if (node(i)%rep>aaa+aa*(icr).and.aa.ne.0) node(i)%rep=aaa+aa*(icr)
        if (node(i)%rep<aaa-aa*(icr).and.aa.ne.0) node(i)%rep=aaa-aa*(icr)
     !   if (aa==0.and.node(i)%rep>delta*icr) node(i)%rep=delta*icr
     !   if (aa==0.and.node(i)%rep<delta*icr*-1) node(i)%rep=delta*icr*-1
    ! RZ end changes

      if (node(i)%rep<0) node(i)%rep=0.0;end if
      if (npag(10)>0) then;ii=10;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%repcel
      aa=nodeoini(10,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%repcel=nodeo(i)%repcel+a

    ! RZ gradual 30-5-17
        if (node(i)%repcel>aaa+aa*(icr).and.aa.ne.0) node(i)%repcel=aaa+aa*(icr)
        if (node(i)%repcel<aaa-aa*(icr).and.aa.ne.0) node(i)%repcel=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%repcel>delta*icr) node(i)%repcel=delta*icr
    !    if (aa==0.and.node(i)%repcel<delta*icr*-1) node(i)%repcel=delta*icr*-1
    ! RZ end changes

      if (node(i)%repcel<0) node(i)%repcel=0.0;end if
      if (npag(11)>0) then;ii=11;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      !aa=node(i)%tor
       aa=nodeoini(11,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%tor=nodeo(i)%tor+a

    ! RZ gradual 30-5-17
        if (node(i)%tor>aaa+aa*(icr).and.aa.ne.0) node(i)%tor=aaa+aa*(icr)
        if (node(i)%tor<aaa-aa*(icr).and.aa.ne.0) node(i)%tor=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%tor>delta*icr) node(i)%tor=delta*icr
    !    if (aa==0.and.node(i)%tor<delta*icr*-1) node(i)%tor=delta*icr*-1
    ! RZ end changes

      if (node(i)%tor<0) node(i)%tor=0.0;end if
      if (npag(12)>0) then;ii=12;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%stor
      aa=nodeoini(12,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%stor=nodeo(i)%stor+a

    ! RZ gradual 30-5-17
        if (node(i)%stor>aaa+aa*(icr).and.aa.ne.0) node(i)%stor=aaa+aa*(icr)
        if (node(i)%stor<aaa-aa*(icr).and.aa.ne.0) node(i)%stor=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%stor>delta*icr) node(i)%stor=delta*icr
    !    if (aa==0.and.node(i)%stor<delta*icr*-1) node(i)%stor=delta*icr*-1
    ! RZ end changes

      if (node(i)%stor<0) node(i)%stor=0.0;end if
      if (npag(13)>0) then;ii=13;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%reqs
      aa=nodeoini(13,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%reqs=nodeo(i)%reqs+a

    ! RZ gradual 30-5-17
        if (node(i)%reqs>aaa+aa*(delta).and.aa.ne.0) node(i)%reqs=aaa+aa*(delta)
        if (node(i)%reqs<aaa-aa*(delta).and.aa.ne.0) node(i)%reqs=aaa-aa*(delta)
    !    if (aa==0.and.node(i)%reqs>delta*icr) node(i)%reqs=delta*icr
    !    if (aa==0.and.node(i)%reqs<delta*icr*-1) node(i)%reqs=delta*icr*-1
    ! RZ end changes

      if (node(i)%reqs<0) node(i)%reqs=0.0;end if
      if (npag(14)>0) then;ii=14;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%ke
      aa=nodeoini(14,2)  ! previous values RZ 30-5-17

   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
    !    if(a<(-delta*icr)*aa) a=(-delta*icr)*aa


      node(i)%ke=nodeo(i)%ke+a

    ! RZ gradual 30-5-17
        if (node(i)%ke>aaa+aa*(icr).and.aa.ne.0) node(i)%ke=aaa+aa*(icr)
        if (node(i)%ke<aaa-aa*(icr).and.aa.ne.0) node(i)%ke=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%ke>delta*icr) node(i)%ke=delta*icr
    !    if (aa==0.and.node(i)%ke<delta*icr*-1) node(i)%ke=delta*icr*-1
    ! RZ end changes

      if (node(i)%ke<0) node(i)%ke=0.0;end if
      if (npag(15)>0) then;ii=15;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii) ;endif;enddo
 
     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%mo
      aa=nodeoini(15,2)  ! previous values RZ 30-5-17


   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
    !    if(a<(-delta*icr)*aa) a=(-delta*icr)*aa

      node(i)%mo=nodeo(i)%mo+a

    ! RZ gradual 30-5-17
        if (node(i)%mo>aaa+aa*(icr).and.aa.ne.0) node(i)%mo=aaa+aa*(icr)
        if (node(i)%mo<aaa-aa*(icr).and.aa.ne.0) node(i)%mo=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%mo>delta*icr) node(i)%mo=delta*icr
    !    if (aa==0.and.node(i)%mo<delta*icr*-1) node(i)%mo=delta*icr*-1
    ! RZ end changes

      if (node(i)%mo<0) node(i)%mo=0.0
      else; if(node(i)%mo==10) node(i)%mo=1000 !!!!!!
      ! RZ NEW, JUST TO ENSURE
end if
      if (npag(16)>0) then;ii=16;a=0 ; do k=1,npag(ii) ; kk=whonpag(ii,k) ; if (gex(i,kk)>0.0d0) then ; 
      a=a+gex(i,kk)*gen(kk)%wa(ii);endif;enddo

     ! Rz 18-5-17 makes node properties to change gradually?????
       icr=icrr ! RZ 30-5-17
      aaa=node(i)%dmo
      aa=nodeoini(16,2)  ! previous values RZ 30-5-17


   !if(a>(delta*icr)*aa) a=(delta*icr)*aa
   !     if(a<(-delta*icr)*aa) a=(-delta*icr)*aa

     node(i)%dmo=nodeo(i)%dmo+a

    ! RZ gradual 30-5-17
        if (node(i)%dmo>aaa+aa*(icr).and.aa.ne.0) node(i)%dmo=aaa+aa*(icr)
        if (node(i)%dmo<aaa-aa*(icr).and.aa.ne.0) node(i)%dmo=aaa-aa*(icr)
    !    if (aa==0.and.node(i)%dmo>delta*icr) node(i)%dmo=delta*icr
    !    if (aa==0.and.node(i)%dmo<delta*icr*-1) node(i)%dmo=delta*icr*-1
    ! RZ end changes

     if (node(i)%dmo<0) node(i)%dmo=0.0;end if

    end if
  end do
!**************************************************************************gradual*******************************
!****************************************************************** PFH 25-05-17 up
  ! if a node gets 3 times the original size of node 1 we fucking kill the program
  if (mod(getot,100).eq.0) then !>>> Is 4-2-14
    if (ffu(7)==1) then  ! if a node gets twice its original size we kill the program ! >>> Is 4-2-14
      do i=1,nd                                                                               ! >>> Is 4-2-14 
        if (node(i)%da>ramax) then                        ! >>> Is 4-2-14
          node(i)%da=ramax
        end if                                                                                ! >>> Is 4-2-14
        !if(node(i)%reqs>1.5)then;print*,"reqs over!",node(i)%reqs;node(i)%reqs=1.5;endif!reqsmax)stop
        !if (node(i)%req>ramax) then   !maximal req now is controlled above in the same sub. !>>Miquel7-8-14
        !  node(i)%req=ramax
        !end if
      end do                                                                                  ! >>> Is 4-2-14
    end if                                                                                    ! >>> Is 4-2-14
  end if                                                                                      ! >>> Is 4-2-14




! Check cpu_time
call cpu_time(cputime)
if((cputime-start_time)>36000)then
  whichend="Over 10 hours 2"
  goto 171 
end if

return !X! THIS IS THE END



! We check in here if all cells are already differentiated, if they are we stop the simulations

!>>> Is 4-4-14

  
171  print *,""
  print *," THIS IS THE END all cells are differentiated and then the simulations stop",trim(carg)//trim(nofi)!trim(noff),
  print *,""
  print*,"reason to end:",whichend
  print*,"total_gabriel",total_gabriel
  call writesnap  
  nofifit=trim(carg)//trim(nofi)
  print*,"nofifit ",nofifit
  call cpu_time(current_time1)
 ! notwin=nofifit(101:103) !filter twins out
 !  notwin=nofifit(60:62)
 ! print*,"notwin",notwin
  !if(notwin=="_10")then
 !   call fit(0,fitelli,nofifit)
 !   open(202,file=trim(carg)//"fitness") 
   ! open(203,file=trim(carg)//"fitnessGene")
   ! open(204,file=trim(carg)//"fitnessAngle")
 !   write(202,*)fitelli    
    !write(203,*)refitGeneN ; write(204,*)refitAngleN
 !   close(202)!;close(203);close(204)   
 ! endif
 ! call cpu_time(current_time2)
 ! print*,"time spend in fitness",nofifit,"-time-",current_time2-current_time1                                                                      !>>> Is 25-2-14
  call cpu_time(stop_time)
  print*,"Time",stop_time-start_time!," STOP by differentation"
  if (len_trim(carg)/=0) then
    print *,trim(carg)!//trim(noff),"ki"
    print *,trim(carg)//trim(nofi),"kii"
    open(23,file=trim(carg)//"t",iostat=i)
    print *,"making...",trim(carg)//"t"
    write(23,*,ERR=46) trim(carg)//trim(nofi)
    flush(23)
    !ret = fsync(fnum(23))
          
    ! Handle possible error
    !if (ret /= 0) stop "Error calling FSYNC"
    open(23,file=trim(carg)//"t",iostat=ii)
    print *,trim(carg)//"t",ii,"iostat"
    read(23,*,END=45,ERR=48) cx
print *,""
print *,cx,"cx HERE"
print *,""
    !close(23)
    call flush(23)
    print *,"done with iostat",ii
    cx="ls -alrt "//trim(carg)//"t" 
    call system(cx) 
  end if
  stop
45 print *,"end of file error"
  stop
46 print *,"error in writing",trim(carg)//"t"
  stop 
48 print *,"other end"
  stop

!>>> Is 4-4-14

end subroutine

!*******************************************************************

subroutine cellbreak   !>>> 17-1-14
  real*8  ax,ay,az,ida
  integer doapo(nd)
  integer nocon(nd,nd)
  integer i,j,k,ii,jj,kk,iii,jjj,kkk,ik

  ! now we check that the cell is not split in two or more parts
  doapo=0
  nocon=0

  do i=1,ncels

    kkk=0
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%marge==0) then ; kkk=ii ; ida=node(kkk)%da ; exit ; end if
    end do

    if (kkk==0) then !it means that the cell has no nucleus and then IT MUST DIE!!!!   !>>> Is 11-6-14
      do j=1,cels(i)%nunodes  !>>> Is 11-6-14
        ii=cels(i)%node(j)    !>>> Is 29-6-14
!print *,ii,i,j,cels(i)%node(:cels(i)%nunodes),ii,"ii",nd
        doapo(ii)=1           !>>> Is 11-6-14
      end do                  !>>> Is 11-6-14
      kkk=1  ! >>> Is 11-6-14
    else
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
        if (sqrt((ax-node(kkk)%x)**2+(ay-node(kkk)%y)**2+(az-node(ii)%z)**2)<node(ii)%da+ida) then
          doapo(ii)=1
        end if
      end do
    end if   ! >>> Is 11-6-14

    doapo(kkk)=1
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      ax=node(ii)%x ; ay=node(ii)%y ; az=node(ii)%z
      ida=node(ii)%da
      do jj=1,cels(i)%nunodes
        if (j==jj) cycle
        iii=cels(i)%node(jj)
        if (sqrt((ax-node(iii)%x)**2+(ay-node(iii)%y)**2+(az-node(iii)%z)**2)<node(iii)%da+ida) then
          nocon(ii,iii)=1
          nocon(iii,ii)=1
          if (doapo(ii)==1) then
            doapo(iii)=1
          else
            if (doapo(iii)==1) doapo(ii)=1
          end if
        end if
      end do
    end do

    do k=1,cels(i)%nunodes/2
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        do jj=1,cels(i)%nunodes
          if (j==jj) cycle
          iii=cels(i)%node(jj)
          if (nocon(ii,iii)==1) then
            if (doapo(ii)==1) then
              doapo(iii)=1
            else
              if (doapo(iii)==1) doapo(ii)=1
            end if
          end if
        end do
      end do
    end do
  end do

  ik=1
  do while(ik<=nd)
    if (doapo(ik)==0) then
      if (node(ik)%tipus>2) then  !we only delete an epithelial node if its altre is also lonely
        call apoptosis(ik)
      else
        if (doapo(node(ik)%altre)==0) then
          call apoptosis(ik)
        end if
      end if
      ik=ik+1
    else
      ik=ik+1   ! I know, it is kind of funky but it should be this way, a loop with nd wont do because nd decreases because of apoptosis
    end if
  end do


end subroutine

!********************************************************************

subroutine polarization
integer:: celd,nnod,tipi,ggr,ccen
real*8::a,b,c,d,e,ax,ay,az,bx,by,bz,cx,cy,cz,ix,iy,iz,alfa,s
      do celd=1,ncels
        tipi=cels(celd)%ctipus
        nnod=cels(celd)%nunodes        
        if (nnod==0) cycle      ! >>> Is 10-5-14
        iy=1d10 ; cx=0d0 ; cy=0d0 ; cz=0d0
	a=cels(celd)%cex ; b=cels(celd)%cey ; c=cels(celd)%cez   

        do i=1,nnod                                                     ! (gen) in the centroid (in the closest node)
          j=cels(celd)%node(i)
          if(node(j)%tipus==1.or.node(j)%tipus==3)then !in epithelial cells, polarity is planar, so we only take one layer of nodes >>>Miquel22-10-13
            d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
          end if
          if(d.le.iy)then;iy=d;ccen=j;endif             
        end do   

        alfa=0.0d0                                                 ! concentration in the central node
        do k=1,npag(nparam_per_node+8)    
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ccen,kk)>0.0d0) then
            alfa=alfa+gex(ccen,kk)*gen(kk)%wa(nparam_per_node+8)   ! wa in units of probability such that it makes things to go from 0 to 1
          end if
        end do  

        ix=0d0 ; iy=0d0 ; iz=0d0                                        ! vector of the gradient within a cell
        do i=1,nnod                                                     
            j=cels(celd)%node(i)
            if(node(j)%tipus==1.or.node(j)%tipus==3)then          
              d=sqrt((node(j)%x-a)**2+(node(j)%y-b)**2+(node(j)%z-c)**2)
              if (d<epsilod) cycle
              d=1d0/d                                                   ! module of radial vectors to get unitary vectors     
              s=0.0d0
              do k=1,npag(nparam_per_node+8)
                kk=whonpag(nparam_per_node+8,k)
                if (gex(j,kk)>0.0d0) then
                  s=s+gex(j,kk)*gen(kk)%wa(nparam_per_node+8)
                end if
              end do
              ix=ix+((node(j)%x-a)*d)*(s-alfa)                   ! and ignore shape/size effects
              iy=iy+((node(j)%y-b)*d)*(s-alfa)
              iz=iz+((node(j)%z-c)*d)*(s-alfa)
            end if
        end do

        if((ix.eq.0).and.(iy.eq.0).and.(iz.eq.0))then            ! if the gene has uniform expresion, the vector is random ! >>>Miguel1-7-14
          call random_number(a)                                  ! >>>Miguel1-7-14
          k=int(a*nvaloq)+1                                      ! >>>Miguel1-7-14
          cels(celd)%polx=particions_esfera(k,1)                 ! >>>Miguel1-7-14
          cels(celd)%poly=particions_esfera(k,2)                 ! >>>Miguel1-7-14
          cels(celd)%polz=particions_esfera(k,3)                 ! >>>Miguel1-7-14
        else                                                     ! >>>Miguel1-7-14
          a=ix**2+iy**2+iz**2 
          if(a==0)then
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0	! unitary resultant vector (gradient polarization)
          else
            d=1d0/sqrt(a)
            cels(celd)%polx=ix*d ; cels(celd)%poly=iy*d ; cels(celd)%polz=iz*d	! unitary resultant vector (gradient polarization)
          end if
          if((ix.eq.0d0).and.(iy.eq.0d0).and.(iz.eq.0d0))then                     ! miguel27-11-13
            cels(celd)%polx=0d0 ; cels(celd)%poly=0d0 ; cels(celd)%polz=0d0
          endif   ! miguel27-11-13
        endif                                                    ! >>>Miguel1-7-14
      end do
end subroutine

!*******************************************************************

subroutine polarizationisaac
integer i,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d,aa,bb,cc 

do i=1,ncels
  a=cels(i)%cex ; b=cels(i)%cey ; c=cels(i)%cez
  sx=0.0d0      ; sy=0.0d0      ; sz=0.0d0
  if (node(cels(i)%node(1))%tipus<3) then
    do j=1,cels(i)%nunodes
      ii=cels(i)%node(j)
      if (node(ii)%tipus==1) then
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+8)*delta
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          d=s/sqrt(aa**2+bb**2+cc**2)
          sx=sx+d*aa ; sy=sy+d*bb ; sz=sz+d*cc
        end if
      end if
    end do
    aa=sqrt(sx**2+sy**2+sz**2)
    if (aa>0.0d0) then
      d=1d0/aa
      cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
    else
      cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
    end if
  else
    if (node(cels(i)%node(1))%tipus==3) then
      do j=1,cels(i)%nunodes
        ii=cels(i)%node(j)
        s=0.0d0
        do k=1,npag(nparam_per_node+8)
          kk=whonpag(nparam_per_node+8,k)
          if (gex(ii,kk)>0.0d0) then
            s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+8)
          end if
        end do
        if (s/=0.0d0) then
          aa=node(ii)%x-a ; bb=node(ii)%y-b ; cc=node(ii)%z-c
          sx=sx+s*aa ; sy=sy+s*bb ; sz=sz+s*cc
        end if
      end do
      aa=sqrt(sx**2+sy**2+sz**2)
      if (aa>0.0d0) then  !epsilod) then
        d=1d0/aa
        cels(i)%polx=sx*d ; cels(i)%poly=sy*d ; cels(i)%polz=sz*d  
      else
        cels(i)%polx=0.0d0 ; cels(i)%poly=0.0d0 ; cels(i)%polz=0.0d0   
      end if
    end if
  end if
end do
end subroutine

!*************************************************************************************

subroutine change_minsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+10)
      kk=whonpag(nparam_per_node+10,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+10)  !wa in units of number of nodes but it can be roughly understood as space-req units
      end if                                           ! THIS WA CAN BE NEGATIVE
    end do
  end do
  s=s*delta
!print *,delta,"ACHTUNG"
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%minsize_for_div=cels(ick)%minsize_for_div+s  !this can make a SUDDEN CHANGE
    if (cels(ick)%minsize_for_div<1) cels(ick)%minsize_for_div=1
!  end if
end do

end subroutine

!*************************************************************************************

!>>> Is 5-2-14
subroutine change_maxsize_for_div  ! updates the size required for dividint according to gene expression nparam_per_node+10
integer ick,j,k,ii,kk
real*8 a,b,c,s,sx,sy,sz,d 

do ick=1,ncels
  s=0.0d0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+15)
      kk=whonpag(nparam_per_node+15,k)
      if (gex(ii,kk)>0.0d0) then
        s=s+gex(ii,kk)*gen(kk)%wa(nparam_per_node+15)  ! THIS WA MAY BE NEGATIVE
      end if
    end do
  end do
  s=s*delta
!  if (s>0.0d0) then
    s=s/cels(ick)%nunodes !this way %minsize is independent of cell size !>>>>Miquel2-12-13
    cels(ick)%maxsize_for_div=cels(ick)%maxsize_for_div+s  !wa in units of space-req roughly this can make a SUDDEN CHANGE
    if (cels(ick)%maxsize_for_div<1) cels(ick)%maxsize_for_div=1
!  end if
end do

!>>> Is 5-2-14

end subroutine


!**************************************************************************************

subroutine emt !epithelial-mensenchymal transitions

               ! this simply makes the tipus equal to 3 but it needs to approach both sides of the originally epithelial cell so that 
               ! the two parts do not drift away since they are more far away than da (they are at reqs normally)
               ! this transition has to be sudden otherwise other existing forces may impede the nodes from getting close enough to each other
               ! that wouldnt be realistic since it would be as if the cell would explode.

integer ick,j,k,ii,kk,iii,jjj,kkk,iv
real*8 a,b,c,d,aa,bb,cc,dd,s

e: do ick=1,ncels
  if (cels(ick)%ctipus==3) cycle
  s=0.0
  do j=1,cels(ick)%nunodes
    ii=cels(ick)%node(j)
    do k=1,npag(nparam_per_node+13)
      kk=whonpag(nparam_per_node+13,k)
      a=gex(ii,kk) !don't change this, weird stuff will happen if you do !>>Miquel18-12-14
      s=s+a*gen(kk)%wa(nparam_per_node+13)
   !if(gex(ii,kk)>epsilod) print*,"emt gene",kk,"node",ii,"cell",ick,"conc",gex(ii,kk)
    end do
  end do
  cels(ick)%temt=cels(ick)%temt+s*delta/real(cels(ick)%nunodes)
  
  if (cels(ick)%temt>=1.0d0) then       ! wa in units of probability: this is an arbitrary value but it is fine since the rest can be re-scaled 
    !so this cell goes to emt
 !print*,"EMT HAPPENS CELL",ick
    cels(ick)%ctipus=3  ! >>> Is 4-1-14
    do jjj=1,cels(ick)%nunodes
      iii=cels(ick)%node(jjj)          
      if (node(iii)%tipus==1) then
        iv=node(iii)%altre
        a=node(iv)%x-node(iii)%x ; b=node(iv)%y-node(iii)%y ; c=node(iv)%z-node(iii)%z 
        d=1d0/sqrt(a**2+b**2+c**2)
        a=a*d ; b=b*d ; c=c*d
        aa=0.5d0*(node(iv)%x+node(iii)%x) ; bb=0.5d0*(node(iv)%y+node(iii)%y) ; cc=0.5d0*(node(iv)%z+node(iii)%z) 
        dd=(node(iii)%da+node(iv)%da)*0.1d0
        node(iii)%x=aa-a*dd 
        node(iii)%y=bb-b*dd 
        node(iii)%z=cc-c*dd
        node(iv)%x=aa+a*dd 
        node(iv)%y=bb+b*dd 
        node(iv)%z=cc+c*dd
        node(iii)%altre=0
        node(iv)%altre=0
        node(iii)%tipus=3
        node(iv)%tipus=3
        !if only one node per cell both nodes of the epithelial cell get a nucleus
        if (ffu(1)==1) then ; node(iii)%marge=0 ; node(iv)%marge=0 ; end if ! >>> Is 10-10-14
      end if
    end do
  end if
end do e

end subroutine

!*******************************************************************************************************

subroutine diffusion_of_reqcr
  real*8 hreqcr(nd),hreqp(nd),hreqc(nd)
  integer i,j,k,ii,jj,kk
  real*8 a,b,c,d

  do i=1,nd
    a=0.0d0 ; b=0.0 ; c=0.0d0
    do ii=1,nneigh(i)
      k=neigh(i,ii)
      if (node(i)%icel/=node(k)%icel) cycle    ! only within the same cell
      if (node(i)%tipus/=node(k)%tipus) cycle  ! only within the same side of the cell
      if (node(i)%tipus>2) cycle               ! only for epithelial cells
      d=dneigh(i,ii)
      a=a+(node(k)%reqcr-node(i)%reqcr) !/(d+1d0)
      !b=b+(node(k)%reqc-node(i)%reqc) !/(d+1d0)
      !c=c+(node(k)%reqp-node(i)%reqp) !/(d+1d0)
    end do
    hreqcr(i)=a*dif_req ! >>> 11-6-14
    !hreqc(i)=b*dif_req  ! >>> 11-6-14
    !hreqp(i)=c*dif_req  ! >>> 11-6-14
  end do  
  do i=1,nd
    node(i)%reqcr=node(i)%reqcr+delta*hreqcr(i)
    !node(i)%reqc=node(i)%reqc+delta*hreqc(i)
    !node(i)%reqp=node(i)%reqp+delta*hreqp(i)
  end do
end subroutine

end module nexus
