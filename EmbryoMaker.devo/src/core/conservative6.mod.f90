module conservative6
! The hereafter following subroutines are called within 8 nested procrustes loops. They are executed for all step sizes that we wish to evaluate. The descendent and ancestorfiles have been read in.

!gfortran -w -fexceptions -fno-underscoring -fbounds-check compare_anton.mod.f90 compare_direct.f90 -o ca.e -lGL -lGLU -lglut
use io
use general
use genetic

integer :: step,boxcount,uz,steps
real*8 :: sizediff ! >>> Is 2 is better 
real*8, allocatable ::boxes_out(:,:)
!real*8 :: sc !,sc0 ! to save scores
real*8 :: precise1,precise2,alpha,beta,gamm !P
real*8 :: maxvo,maxv,maxxc,maxyc,maxzc,minxc,minyc,minzc,innd !P
real*8, allocatable :: c0x(:),c0y(:),c0z(:)
integer, allocatable :: boxes0(:,:,:,:),boxcount0(:)
real*8,allocatable :: points0(:,:)!,points(:,:)
real*8 :: oo,pp ! extension old
integer :: nd0,boxnum,nddd,ndtt,pasadas,ndop
integer :: mxc,myc,mzc,nx,ny,nz !P
integer :: boxn=300,boxm=300,boxs=3 ! maximal and minimal number of partitioning steps
real*8 :: minprecise,maxprecise,smaller ! minimal degree of precision
! character*140 pernofil,pernofila
 character*140 ancestor,descens ! ancestor and descendent files
!logical :: ex
real*8 ::  x,y,zc,aao,sizedifi
real*8,dimension(3) :: medi
character*140 ancestori,ancestorii
character*1 jabato
contains

subroutine conservative(sc,pernofil,pernofil2)!(indiv,sc,pernofil,pernofila,gd1) !(sc,pernofil)   ! conservative(sc)

character*140 pernofil,pernofil2

real*8 :: sc,ges,o
real*8,allocatable :: get(:)
integer :: indiv,u
real*8 :: gd1,gd2


descens=trim(pernofil)
ancestor=trim(pernofil2)


print*, "THESE ARE THE FILES:", descens, ancestor

!call boxnumbers(ancestor,descens,boxnum)
!open(39,file=trim(descens)//"_EMD")
!print*, boxnum
smaller=0.7
pasadas=0
if(allocated(c0x))then
  deallocate(c0x)
  deallocate(c0y)
  deallocate(c0z)
end if
allocate(c0x(1),c0y(1),c0z(1))
c0x=0d0; c0y=0d0; c0z=0d0

!call ancestor_define_boxes(ancestor,oo,int(boxnum*smaller),sizediff) ! this should be called only once during reva

!++++++++++
! starts here!!
print*, "Now read in the new morphology!" ! strangest error ever: if I do not print anything here, there will be an input error

nd=0
ndtt=0
call iniread
call readsnap(ancestor)
print*,"nd&&",nd
do j=1,nd
 if(node(j)%tipus==2) then
  ndtt=ndtt+1
 endif
end do
ndop=nd
109 print*, "READIN DONE", ancestor, " POINTS: ", nd, ndtt
if(allocated(points0)) deallocate(points0)
allocate(points0(nd,3))

do j=1,nd
 if(node(j)%tipus==2) then
  !ndtt=ndtt+1
  points0(j,1)=node(j)%x
  points0(j,2)=node(j)%y
  points0(j,3)=node(j)%z
 else
  points0(j,:)=4004
 endif
end do

!PFH print
! open(311,file="ini_optimum.csv")
!   do j=1,ndop
!    if(points0(j,1)>4000)cycle
!    write(311,*) points0(j,1),",",points0(j,2),",",points0(j,3)
!   enddo
! close(311)

!nd=nd0
!ndtt=nd
print*, "NDTT", ndtt
steps=1

 minxc=2000;minyc=2000;minzc=2000;maxxc=-1*huge(maxxc);maxyc=-1*huge(maxyc);maxzc=-1*huge(maxzc) !huge(minxc) huge(minyc)  huge(minzc)
 maxvo=0.0d0
 c0x=0.0d0; c0y=0.0d0; c0z=0.0d0

 do j=1,ndop ! save the values of the minimal nodes of x,y,z
   if(points0(j,1)>4000)cycle
   if(points0(j,1).lt.minxc) minxc=points0(j,1)
   if(points0(j,2).lt.minyc) minyc=points0(j,2)
   if(points0(j,3).lt.minzc) minzc=points0(j,3)
   if(points0(j,1).gt.maxxc) maxxc=points0(j,1)
   if(points0(j,2).gt.maxyc) maxyc=points0(j,2)
   if(points0(j,3).gt.maxzc) maxzc=points0(j,3)
 end do 

maxvo=maxxc-minxc ! now search for the maximal extension in whatever direction; it is important for the absolute size
 if(maxvo.lt.maxyc-minyc) maxvo=maxyc-minyc
 if(maxvo.lt.maxzc-minzc) maxvo=maxzc-minzc
print*,"nd",nd,"ndtt",ndtt
 aao=0d0; pp=0d0
!print*,"pp1",pp
 do nodj=1,ndop
   if(points0(nodj,1)>4000)cycle
   c0x=c0x+(points0(nodj,1))-minxc; c0y=c0y+(points0(nodj,2))-minyc; c0z=c0z+(points0(nodj,3))-minzc  ! >>> RZ 29-06 centroid 
   do j=1,ndop
    if(points0(j,1)>4000)cycle
    if(nodj==j)cycle
     p=(points0(nodj,1)-points0(j,1))**2+(points0(nodj,2)-points0(j,2))**2+(points0(nodj,3)-points0(j,3))**2!;print*,"p**",p
     if(p>pp) pp=p
   end do
 end do
!print*,"pp2",pp
 pp=sqrt(pp)
!print*,"pp3",pp
 sizedifi=pp

 c0x=c0x/real(ndtt); c0y=c0y/real(ndtt); c0z=c0z/real(ndtt) ! >>> RZ 29-06 centroid 

 !steps=boxnum2  
 !  boxcount0(steps)=sum(boxes0(:,:,:,steps)) ! the total number of occupied boxes
   c0x(steps)=c0x(steps)/maxvo*steps; c0y(steps)=c0y(steps)/maxvo*steps; c0z(steps)=c0z(steps)/maxvo*steps ! >> RZ 29-06 centroid

do nodj=1,ndop
  if(points0(nodj,1)>4000)cycle
  points0(nodj,1)=(points0(nodj,1)-minxc)/real(maxvo)
  points0(nodj,2)=(points0(nodj,2)-minyc)/real(maxvo)
  points0(nodj,3)=(points0(nodj,3)-minzc)/real(maxvo)
end do

call iniread
call readsnap(descens)

nddd=0
print*,"pp4",pp
do j=1,nd
 if(node(j)%tipus==2) then
  nddd=nddd+1
 endif
end do

print*, "READIN2 DONE", descens, " POINTS: ", nd, nddd

call conservative_main(sc,descens,oo,int(boxnum*smaller),sizediff) 
!PFH print
! open(35,file="line179_final_optimum.csv")
!   do j=1,ndop
!    if(points0(j,1)>4000)cycle
!    write(35,*) points0(j,1),",",points0(j,2),",",points0(j,3)
!   enddo
! close(35)
!print*,"ver final_optimum y otros"
!235 print*,"continue? y or n"
!read(*,*)jabato
!if(jabato=="n")goto 235 
 

229 print*, "THIS IS YOUR DISTANCE:", sc

end subroutine conservative
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine conservative_main(sc,descensi,ooo,boxnum1,sizedif)

integer :: i,rx,ry,dx,dy,dz,sizes,runthrough,boxnum1,preprecise
integer :: mxc,myc,mzc,nx,ny,nz,o1,o2,o3,ka,ke,uy
real*8 :: sizedif,q,qq,ku
real*8 :: sco,scn,scm,sc00,alpha,beta,sizescale,pi,sc,cx,cy,cz,ddx,ddy,ddz,precise2,precise3,nn,ooo,stepminp=0.000001
real*8, dimension(boxn) :: score,scor
real*8, dimension(6,3,3) :: keep
real*8, dimension(nd,3) :: des,des1,dess,des_a,des_b,ori_des,points
real*8 :: stepmin=0.0001 ! minimal score difference to allow for further partitioning
integer, allocatable :: boxes1(:,:,:)
real*8  px,py,pz,ppx,ppy,ppz,ao,bo,co,saves
integer,dimension(3,2) :: boundzz
 character*140 descensi

real*8 :: zzz,zz0,zz1

ku=0.5d0!0.5
print*, nddd, "########"

! NEW
step=1; boxnum1=1

precise3=0.25
pi=3.14159265359
minprecise=0.005;maxprecise=1E-8
des=0.0d0; dess=0.0d0; ori_des=0.0d0; des1=0.0d0
precise2=1; preprecise=12; score=huge(score)
step=boxs ! START BY 3
innd=1.0d0/nddd ! >> RZ 29-06 runs faster

! THESE VALUES CAN BE MODIFIED AD LIBITUM
runthrough=1 !if 0, it only uses the "small procrustes loop", overlay plus rotation (no further translations)
precise2=1 ! starting precision
preprecise=12 ! maximal precision of the rotations in the "small procrustes loop"
!boxm and boxn to be set above: number of box divisions

call iniread
call readsnap(descensi)
!nddd=nd

do j=1,nd
  if(node(j)%tipus==2) then
!  nddd=nddd+1
   points(j,1)=node(j)%x
   points(j,2)=node(j)%y
   points(j,3)=node(j)%z
  else
   points(j,:)=4004
 endif
end do
!PFH print
! open(309,file=trim(descensi)//"_pre_procrustes.csv")
!   do i=1,nd
!    if(points(i,1)>4000)cycle
!    write(309,*) points(i,1),",",points(i,2),",",points(i,3)
!   enddo
! close(309)

! 2. Identify the matrices' longest axes
a=0d0; b=1E-8; c=1E-8; d=1E-8; o=0d0; qq=0d0
!print*,"**points 234**",points
do i=1,nd
 if(points(i,1)>4000)cycle
  do j=1,nd
    if(i==j) cycle
    if(points(j,1)>4000)cycle
    q=(points(j,1)-points(i,1))**2+(points(j,2)-points(i,2))**2+(points(j,3)-points(i,3))**2
    if(q>qq) qq=q
    o=points(j,1)-points(i,1)
    if(o>a) a=o
    o=points(j,2)-points(i,2)
    if(o>a) a=o
    o=points(j,3)-points(i,3)
    if(o>a) a=o
  end do
  o=points(i,1)
  if(o<b) b=o
  o=points(i,2)
  if(o<c) c=o
  o=points(i,3)
  if(o<d) d=o
end do
o=a
q=sqrt(qq)

do i=1,nd
 if(points(i,1)>4000)then;ori_des(i,:)=4004;cycle;endif
 ori_des(i,1)=points(i,1)-b; ori_des(i,2)=points(i,2)-c; ori_des(i,3)=points(i,3)-d 
enddo
!print*,"ori_des 263",ori_des!
!	open(15, file="FF")
!     		do uy=1,nd
!                        if(ori_des(uy,1)>4000)cycle
!     			write(15,*) ori_des(uy,1:3)
!     		end do
!     	close(15)

! 1.Partitions loop

step=1 

 389    sc=huge(sc)*0.1;sco=huge(sco)*0.3

    des1(:,1)=(ori_des(:,1))
    des1(:,2)=(ori_des(:,2))
    des1(:,3)=(ori_des(:,3))
!print*,"des1 279",des1
sco=0d0!score(step)  ! >> RZ 29-06 new
des=des1
sc=0d0

! for very similar or identical object that are not turned (the "usual" case), it makes no sense to go through the procrustes loops

!print*, "PPP", pp, q

  do i=1,nd
    if(des(i,1)>4000)cycle
    do j=1,3
      des(i,j)=des(i,j)/q!*(pp/q)
    end do
  end do

!       open(19, file="line296")
!     		do uy=1, nd
!                        if(des(uy,1)>4000)cycle
!     			write(19,*) des(uy,1:3)
!     		end do
!     	close(19)

  des1=des
!print*,"des 296",des
  call find_centroid(des1,cx,cy,cz) ! >>> RZ 29-06 centroid 

  c0x(step)=c0x(step)-cx; c0y(step)=c0y(step)-cy; c0z(step)=c0z(step)-cz !>>> RZ 29-06
 
  ddx=c0x(step); ddy=c0y(step); ddz=c0z(step)

!open(16,file="CENTER")
!write(16,*) ddx, ddy, ddz, cx, cy, cz
close(16)


des=des1
!print*,"des 309",des
do i=1,nd
 if(des(i,1)>4000)cycle
 des(i,1)=des(i,1)+ddx; des(i,2)=des(i,2)+ddy; des(i,3)=des(i,3)+ddz
end do
!open(291, file="line322")
!     		do uy=1, nd
!                        if(des(uy,1)>4000)cycle
!     			write(291,*) des(uy,1:3)
!     		end do
!     	close(291)
des1=des
call compare(des,sc)
!print*,"sc***337",sc
!print*,"des 315",des1
call find_centroid(des1,cx,cy,cz)
!open(299, file="line338")
!     		do uy=1, nd
!                        if(des(uy,1)>4000)cycle
!     			write(299,*) des(uy,1:3)
!     		end do
!     	close(299)
!nddd=nd
  ! pre-rotate

  call compare(des,sc)
!print*,"sc***348",sc
  score(step)=sc
  scor(step)=sc
  scm=sc
  des1=des

     	!	open(15, file="DD")
     	!	do uy=1,nddd
     	!		write(15,*) des1(uy,1:3)
     	!	end do
     	!	close(15)


  gamm=0d0
  if(ku==1)then
    ka=-7; ke=8
  elseif(ku==0.5)then
    ka=-15;ke=16
  endif

  
  keep=0d0
  do i=0,preprecise!*3 ! how precise in the beginning
  precise2=1d0/(real(i*0.5+2))
  alpha=0d0; beta=0d0; gamm=0d0
  do rx=ka,ke   ! rotation in x-direction
    alpha=rx*0.25*pi*precise2*ku
    do ry=ka,ke   ! rotation in y-direction
      beta=ry*0.25*pi*precise2*ku
     ! do rz=ka,ke   ! rotation in z-direction
       ! gamm=rz*0.25*pi*precise2*ku
        des=des1
        !call rotater3(des,alpha,beta,gamm)
        call rotater(des,alpha,beta)
        !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
        call compare(des,sc)
        if(sc.le.scm)then
          scm=sc
          keep(4,1,1:3)=alpha; keep(5,1,1:3)=beta; keep(6,1,1:3)=gamm
          o2=3
        end if
     ! end do
    end do
  end do
  !des1=des
  !call rotater3(des1,keep(4,1,1),keep(5,1,1),keep(6,1,1))

     !           open(13, file="DD1",position="append")
     !		do uy=1,nddd
     !			write(13,*) des1(uy,1:3), i
     !		end do
     !		close(13)

call rotater(des1,keep(4,1,1),keep(5,1,1))

  des=des1
 ! des1=des

     	!	open(14, file="DD2",position="append")
     	!	do uy=1,nddd
     	!		write(14,*) des1(uy,1:3), i
     	!	end do
     	!	close(14)

  end do

dd=0d0
     !		open(15, file="CC")
     		do uy=1, nd
                       if(des1(uy,1)>4000)cycle
     	!		write(15,*) des1(uy,1:3)
                       do uz=1,nd
                         if(uz==uy) cycle
                         if(des1(uz,1)>4000)cycle
                         d=(des1(uy,1)-des1(uz,1))**2+(des1(uy,2)-des1(uz,2))**2+(des1(uy,3)-des1(uz,3))**2
                         if(d>dd) dd=d
                       end do
     		end do
     	!	close(15)

!print*, "DC", sqrt(dd)
!stop

    score(step)=scm
    scn=scm
    saves=scm
if(runthrough==0) goto 129
  sc=0d0
  alpha=0d0; beta=0d0; gamm=0d0
  precise2=precise2*2.0 ! we do not need to start from the original precise any more
  ! 2.Precision loop/Initial conditions
    do while(((precise2.gt.minprecise).or.(abs(sco-sc)/(sco+sc).gt.stepminp)).and.(precise2.gt.maxprecise))  ! precision !
     sco=score(step)  ! >> RZ 29-06 new
    ! Initial conditions: 1.No change; translate scaled with centroids, 2.Overlay centroids, scale absolute???
   ! 3. PROCRUSTES subroutines call 
      !. Translation first
     sc00=sc
     keep=0d0
     keep(1:3,:,:)=0d0; keep(4:5,:,:)=0d0; keep(6,:,:)=1d0
     scn=huge(scn); scm=huge(scm)
     do dx=-2,2   ! translation in x-direction
       ddx=dx*0.5*precise2
       do dy=-2,2   ! translation in y-direction
         ddy=dy*0.5*precise2
         do dz=-2,2   ! translation in y-direction
           ddz=dz*0.5*precise2
           des=des1
           call translate(des,ddx,ddy,ddz)
           !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
           call compare(des,sc)
           if(sc.le.scm)then
             scm=sc
             keep(1,1:2,1)=ddx; keep(2,1:2,1)=ddy; keep(3,1:2,1)=ddz
             o1=1
           end if
         end do
       end do
     end do
     scn=scm  !!!!!!!!!!!
     ! T - R - S
     des=des1
     call translate(des,keep(1,1,1),keep(2,1,1),keep(3,1,1))
     dess=des
     do rx=-3,4   ! rotation in x-direction
       alpha=rx*0.5*pi*precise2
       do ry=-3,4   ! rotation in y-direction
         beta=ry*0.5*pi*precise2
         des=dess
         call rotater(des,alpha,beta)
         !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
         call compare(des,sc)
         if(sc.le.scm)then
           scm=sc
           keep(4,1,1)=alpha; keep(5,1,1)=beta
           o2=3
         end if
       end do
     end do
     des=dess
     call rotater(des,keep(4,1,1),keep(5,1,1))
     dess=des
     do sizes=-2,2
       sizescale=1.0d0+sizes*0.25*precise2 !>>> RZ 29-06 moved here 
       des=dess
       call scaling(des,sizescale)
       !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
       call compare(des,sc)
       if(sc.le.scm)then
         scm=sc
         keep(6,1,1)=sizescale
         o3=2
       end if
     end do
     des=dess
     call scaling(des,keep(6,1,1))

     ! T - S - R
     des=des1
     call translate(des,keep(1,2,1),keep(2,2,1),keep(3,2,1))
!goto 335
     dess=des
     do sizes=-2,2
       sizescale=1.0d0+sizes*0.25*precise2 !>>> RZ 29-06 moved here 
       des=dess
       call scaling(des,sizescale)
       !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
       call compare(des,sc)
       if(sc.le.scn)then
         scn=sc
         keep(6,2,1)=sizescale
         o1=1; o2=2
       end if
     end do
     des=dess
     call scaling(des,keep(6,2,1))
     dess=des
     do rx=-3,4   ! rotation in x-direction
       alpha=rx*0.5*pi*precise2
       do ry=-3,4   ! rotation in y-direction
         beta=ry*0.5*pi*precise2
         des=dess
         call rotater(des,alpha,beta)
         !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
         call compare(des,sc)
         if(sc.le.scn)then
           scn=sc
           keep(4,2,1)=alpha; keep(5,2,1)=beta
           o2=2;o3=3
         end if
       end do
     end do
     des=dess
     call rotater(des,keep(4,2,1),keep(5,2,1))
     if((scn.le.scm).and.(scn.lt.sc00))then
       score(step)=scn
       keep(:,3,1)=keep(:,2,1)
     elseif((scm.le.scn).and.(scm.lt.sc00))then
       score(step)=scm
       keep(:,3,1)=keep(:,1,1)
     endif
     ! Scaling first
!goto 335
     scn=huge(scn); scm=huge(scm)
     do sizes=-2,2
       sizescale=1.0d0+sizes*0.25*precise2 !>>> RZ 29-06 moved here 
       des=des1
       call scaling(des,sizescale)   
       !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
       call compare(des,sc)
       if(sc.le.scm)then
         scm=sc
         keep(6,1:2,2)=sizescale
         o1=2
       end if
     end do
     ! S - T - R
     des=des1
     call scaling(des,keep(6,1,2))
     dess=des
     do dx=-2,2   ! translation in x-direction
       ddx=dx*0.5*precise2
       do dy=-2,2   ! translation in y-direction
         ddy=dy*0.5*precise2
         do dz=-2,2   ! translation in y-direction
           ddz=dz*0.5*precise2
           des=dess
           call translate(des,ddx,ddy,ddz)
           !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
           call compare(des,sc)
           if(sc.le.scm)then
             scm=sc
             keep(1,1,2)=ddx; keep(2,1,2)=ddy; keep(3,1,2)=ddz
             o2=1; o1=2
           end if
         end do
       end do
     end do
     des=dess
     call translate(des,keep(1,1,2),keep(2,1,2),keep(3,1,2))
     dess=des
     do rx=-3,4   ! rotation in x-direction
       alpha=rx*0.5*pi*precise2
       do ry=-3,4   ! rotation in y-direction
         beta=ry*0.5*pi*precise2
         des=dess
         call rotater(des,alpha,beta)
         !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
         call compare(des,sc)
         if(sc.le.scm)then
           scm=sc
           keep(4,1,2)=alpha; keep(5,1,2)=beta
           o3=3; o2=1
         end if
       end do
     end do
     ! S - R - T
     des=des1
     call scaling(des,keep(6,2,2))
     dess=des
     do rx=-3,4   ! rotation in x-direction
       alpha=rx*0.5*pi*precise2
       do ry=-3,4   ! rotation in y-direction
         beta=ry*0.5*pi*precise2
         des=dess
         call rotater(des,alpha,beta)
         !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
         call compare(des,sc)
         if(sc.le.scn)then
           scn=sc
           keep(4,2,2)=alpha; keep(5,2,2)=beta
           o1=2; o2=3
         end if
       end do
     end do
     des=dess
     call rotater(des,keep(4,2,2),keep(5,2,2))
     dess=des
     do dx=-2,2   ! translation in x-direction
       ddx=dx*0.5*precise2
       do dy=-2,2   ! translation in y-direction
         ddy=dy*0.5*precise2
         do dz=-2,2   ! translation in y-direction
           ddz=dz*0.5*precise2
           des=dess
           call translate(des,ddx,ddy,ddz)
          ! call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
           call compare(des,sc)
           if(sc.le.scn)then
             scn=sc
             keep(1,2,2)=ddx; keep(2,2,2)=ddy; keep(3,2,2)=ddz
             o2=3; o3=1
           end if
         end do
       end do
     end do
     if(scn.le.score(step))then
       score(step)=scn
       keep(:,3,2)=keep(:,2,2)
     else if(scm.le.score(step))then
       score(step)=scm
       keep(:,3,2)=keep(:,1,2)
     endif
     ! Rotation first

     scn=huge(scn); scm=huge(scm)
     do rx=-3,4   ! rotation in x-direction
       alpha=rx*0.5*pi*precise2
       do ry=-3,4   ! rotation in y-direction
         beta=ry*0.5*pi*precise2
         des=des1
         call rotater(des,alpha,beta)
         !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
         call compare(des,sc)
         if(sc.le.scm)then
           scm=sc
           keep(4,1:2,3)=alpha; keep(5,1:2,3)=beta
           o1=3
         end if
       end do
     end do
     ! R - T - S
     des=des1
     call rotater(des,keep(4,1,3),keep(5,1,3))
     dess=des
     do dx=-2,2   ! translation in x-direction
       ddx=dx*0.5*precise2
       do dy=-2,2   ! translation in y-direction
         ddy=dy*0.5*precise2
         do dz=-2,2   ! translation in y-direction
           ddz=dz*0.5*precise2
           des=dess
           call translate(des,ddx,ddy,ddz)
          ! call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
           call compare(des,sc)
           if(sc.le.scm)then
             scm=sc
             keep(1,1,3)=ddx; keep(2,1,3)=ddy; keep(3,1,3)=ddz
             o2=1;o1=3
           end if
         end do
       end do
     end do
     des=dess
     call translate(des,keep(1,1,3),keep(2,1,3),keep(3,1,3))
     dess=des
     do sizes=-2,2
       sizescale=1.0d0+sizes*0.25*precise2 !>>> RZ 29-06 moved here 
       des=dess
       call scaling(des,sizescale)
      ! call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
       call compare(des,sc)
       if(sc.le.scm)then
         scm=sc
         keep(6,1,3)=sizescale
         o3=2;o2=1
       end if
     end do
     ! R - S - T
     des=des1
     call rotater(des,keep(4,2,3),keep(5,2,3))
     dess=des
     do sizes=-2,2
       sizescale=1.0d0+sizes*0.25*precise2 !>>> RZ 29-06 moved here 
       des=dess
       call scaling(des,sizescale)
       !call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
       call compare(des,sc)
       if(sc.le.scn)then
         scn=sc
         keep(6,2,3)=sizescale
         o2=2;o1=3
       end if
     end do
     des=dess
     call scaling(des,keep(6,2,3))
     dess=des
     do dx=-2,2   ! translation in x-direction
       ddx=dx*0.5*precise2
       do dy=-2,2   ! translation in y-direction
         ddy=dy*0.5*precise2
         do dz=-2,2   ! translation in y-direction
           ddz=dz*0.5*precise2
           des=des1
           call translate(des,ddx,ddy,ddz)
          ! call define_boxes(des,boxes1,step,boundzz,precise2)  ! the boxes of the descendant based on the box sizes of the ancestral grid
           call compare(des,sc)
           if(sc.le.scn)then
             scn=sc
             keep(1,2,3)=ddx; keep(2,2,3)=ddy; keep(3,2,3)=ddz
             o3=1;o2=2
           end if
         end do
       end do
     end do
     if(scn.le.score(step))then
       score(step)=scn
       keep(:,3,3)=keep(:,2,3)
     else if(scm.le.score(step))then
       score(step)=scm
       keep(:,3,3)=keep(:,1,3)
     endif
 ! so, save the optimal configuration and recapitulate it
     335 des=des1
     if(saves<score(step))then
       score(step)=saves
       goto 231
     end if

     if(o1==1)then
       call translate(des,keep(1,3,1),keep(2,3,1),keep(3,3,1))
     elseif(o1==3)then
       call rotater(des,keep(4,3,3),keep(5,3,3))
     else
       call scaling(des,keep(6,3,2))
     endif
     if(o2==1)then
       call translate(des,keep(1,3,o1),keep(2,3,o1),keep(3,3,o1))
     elseif(o2==3)then
       call rotater(des,keep(4,3,o1),keep(5,3,o1))
     else
       call scaling(des,keep(6,3,o1))
     endif
     if(o3==1)then
       call translate(des,keep(1,3,o1),keep(2,3,o1),keep(3,3,o1))
     elseif(o3==3)then
       call rotater(des,keep(4,3,o1),keep(5,3,o1))
     else
       call scaling(des,keep(6,3,o1))
     endif
     des1=des
		!if(step==16)then
		!print*, "test", trim(descensi)
  !   		open(15, file="procrustes_res")
  !   		do uy=1, nd
  !   			write(15,*) des1(uy,1:3)
  !   		end do
  !   		close(15)
     		!endif
    231 precise2=precise2*0.5  ! >> RZ 29-06 new

  end do ! >> RZ 29-06 new

  mxc=boundzz(1,1); nx=boundzz(1,2)
  myc=boundzz(2,1); ny=boundzz(2,2)
  mzc=boundzz(3,1); nz=boundzz(3,2)

129    if(score(step).gt.scor(step)) score(step)=scor(step)

   !step=step+1
  ! print*, "NEXT STEP SIZE: ", step, " SCORE: ", score(step-1)

388 sc=0d0
!end do
!print*, score
sc=score(boxnum1)
!if(step-1.lt.boxn) score(step:boxn)=score(step-1)
985 print*, "DONE"!do i=3,boxn; if(score(i)>1) score(i)=0d0; sc=sc+score(i); print*, i, sc, score(i); enddo
!print*,"pre-final sc",sc
call comparefinal(des,sc);print*,"final sc",sc !l compare44
!sc=sc/(boxn-2)
!if(sc<4d-7) sc=0d0
 !sc=sc*1d12
 !pfh print
! open(36,file="line835_final.csv")
!   do i=1,nd
!    if(des(i,1)>4000)cycle
!    !if(des(i,1)<-1)cycle
!    write(36,*) des(i,1),",",des(i,2),",",des(i,3)
!   enddo
! close(36)

! close(39)

302 print*, "THE SCORE", sc

end subroutine conservative_main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine scaling(des,sizescales)

real*8, dimension(nd,3) :: des
real*8 :: sizescales
real*8 :: cx,cy,cz,bb,ss

call find_centroid(des,cx,cy,cz)
do i=1,nd
  if(des(i,1)>4000)then;cycle;else
    des(i,1)=(des(i,1)-cx)*sizescales+cx 
    des(i,2)=(des(i,2)-cy)*sizescales+cy
    des(i,3)=(des(i,3)-cz)*sizescales+cz
  endif
enddo
! BTW: This was wrong, because it should not rescale AND translate at the same time.

end subroutine scaling
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine translate(des,ddx,ddy,ddz)

real*8, dimension(nd,3) :: des
real*8 :: ddx, ddy, ddz !>>> RZ 29-06 REAL instead of INTEGER

! Achtung!: minxyz,maxv have been calculated based on the original coordinates in the define_boxes subroutine. They denote the extreme coordinates of the original descendent morphology.
do i=1,nd
  if(des(i,1)>4000)then;cycle;else
   des(i,1)=des(i,1)+ddx
   des(i,2)=des(i,2)+ddy
   des(i,3)=des(i,3)+ddz
  endif
enddo

end subroutine translate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotater(des3,alpha,beta)

real*8 :: cx,cy,cz
integer i
real*8, dimension(nd,3) :: des3
real*8  px,py,pz,ppx,ppy,ppz
real*8 :: alpha,beta

!print*, "###", nd

!open(11,file="II")
!do i=1,nd
!write(11,*) des3(i,1:3), "0"
!end do

call find_centroid(des3,cx,cy,cz) ! since rotations are around (0,0,0)

!do i=1,nd
!write(11,*) des3(i,1:3), "3"
!end do
!if(2==3)then
!alpha=0.5; beta=0.5
do i=1,nd
  if(des3(i,1)>4000)cycle
  px=des3(i,1)-cx; py=des3(i,2)-cy; pz=des3(i,3)-cz
  ppy=py*cos(alpha)-pz*sin(alpha)
  ppz=pz*cos(alpha)+py*sin(alpha)
  ppx=px*cos(beta)+ppz*sin(beta)
  ppz=ppz*cos(beta)-px*sin(beta)
!write(11,*) des3(i,1:3), "1"
  des3(i,1)=ppx+cx; des3(i,2)=ppy+cy; des3(i,3)=ppz+cz
!!write(11,*) des3(i,1:3), "2"
end do
close(11)
!stop
!end if

end subroutine rotater
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rotater3(des,alpha,beta,gamm)

real*8 :: cx,cy,cz
integer i
real*8, dimension(nd,3) :: des
real*8  px,py,pz,ppx,ppy,ppz
real*8 :: alpha,beta,gamm

call find_centroid(des,cx,cy,cz) ! since rotations are around (0,0,0)

!if(9==0)then
!do i=1,nd
!  px=des(i,1)-cx; py=des(i,2)-cy; pz=des(i,3)-cz
!  ppy=py*cos(alpha)-pz*sin(alpha)
!  ppz=pz*cos(alpha)+py*sin(alpha)
!  ppx=px*cos(beta)+ppz*sin(beta)
!  ppz=ppz*cos(beta)-px*sin(beta)
!  ppx=ppx*cos(gamm)-ppy*sin(gamm)
!  ppy=ppx*sin(gamm)+ppy*cos(gamm)
!  des(i,1)=ppx+cx; des(i,2)=ppy+cy; des(i,3)=ppz+cz
!end do
!endif

do i=1,nd
  if(des(i,1)>4000)cycle
  px=des(i,1)-cx; py=des(i,2)-cy; pz=des(i,3)-cz
  ppy=py*cos(alpha)-pz*sin(alpha)
  ppz=pz*cos(alpha)+py*sin(alpha)

  ppx=px*cos(beta)+ppz*sin(beta)
  ppz=ppz*cos(beta)-px*sin(beta)
  px=ppx; py=ppy
  ppx=px*cos(gamm)-py*sin(gamm)
  ppy=px*sin(gamm)+py*cos(gamm)
  des(i,1)=ppx+cx; des(i,2)=ppy+cy; des(i,3)=ppz+cz
end do

end subroutine rotater3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_centroid(des,cx,cy,cz)

real*8 :: cx,cy,cz
real*8, dimension(nd,3) :: des

 cx=0.0d0; cy=0.0d0; cz=0.0d0

 do i=1,nd  ! calculate the centroid 
   if(des(i,1)>4000)cycle
   cx=cx+des(i,1); cy=cy+des(i,2); cz=cz+des(i,3)
 end do
 cx=cx*innd; cy=cy*innd; cz=cz*innd

end subroutine find_centroid 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine define_boxes(des,boxes1,stepi,boundz,prcs2)

integer :: j,nnax,stepi
integer, dimension(3,2) :: boundz
real*8 :: minx2,miny2,minz2,maxx2,maxy2,maxz2,mmmm,m1,m2,m3,mmax,prcs2
real*8, dimension(nd,3) :: des
integer, allocatable :: boxes1(:,:,:)

minx2=huge(minx2);miny2=huge(miny2);minz2=huge(minz2);maxx2=-1*huge(maxx2);maxy2=-1*huge(maxy2);maxz2=-1*huge(maxz2)
minx2=100000000d0;miny2=100000000d0;minz2=100000000d0
maxx2=-100000000d0;maxy2=-100000000d0;maxz2=-100000000d0

 do j=1,nddd ! save the values of the minimal nodes of x,y,z
   if(des(j,1).lt.minx2) minx2=des(j,1)
   if(des(j,2).lt.miny2) miny2=des(j,2)
   if(des(j,3).lt.minz2) minz2=des(j,3)
   if(des(j,1).gt.maxx2) maxx2=des(j,1)
   if(des(j,2).gt.maxy2) maxy2=des(j,2)
   if(des(j,3).gt.maxz2) maxz2=des(j,3)
 end do

 mxc=int(minx2); nx=int(maxx2+1)
 myc=int(miny2); ny=int(maxy2+1)
 mzc=int(minz2); nz=int(maxz2+1)

 mmax=nx; nnax=1
 m1=maxx2-minx2; m2=maxy2-miny2; m3=maxz2-minz2
 if(m2.gt.m1)then; mmax=ny; nnax=2; endif
 if((m3.gt.m2).and.(m3.gt.m1))then; mmax=nz; nnax=3; endif
 mmmm=tiny(mmmm)

if(allocated(boxes1)) deallocate(boxes1)
allocate(boxes1(mxc:nx,myc:ny,mzc:nz))
boxes1=0

do i=1,nddd !>>> RZ 29-06 This one was wrong

  boxes1(int(des(i,1)),int(des(i,2)),int(des(i,3)))=1 !>>> RZ 29-06

end do !>>> RZ 29-06

557 i=1

boundz(1,1)=mxc; boundz(1,2)=nx
boundz(2,1)=myc; boundz(2,2)=ny
boundz(3,1)=mzc; boundz(3,2)=nz

end subroutine define_boxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compare2(dese,sce)
real*8 :: sce,aaa,a,EMD
!integer :: ndop
real*8 , dimension(nd,3) :: dese
 character*10 pasadach
!integer*8 :: ndtt
sce=0
!a=huge(a)
EMD=0.0
pasadas=pasadas+1
!ndop=ndtt
!write(pasadach,'(I10)')pasadas
  
! open(38,file=trim(descensi)//"_final_"//trim(adjustl(pasadach)))
!   do i=1,nd
!    write(38,*) dese(i,1),dese(i,2),dese(i,3)
!   enddo
! close(38)

  
  do i=1,ndop
   a=1d4
   if(points0(i,1)>4000)cycle
   do j=1,nd!dd
    if(dese(j,1)>4000)cycle
    aaa=(points0(i,1)-dese(j,1))**2+(points0(i,2)-dese(j,2))**2+(points0(i,3)-dese(j,3))**2
    if(aaa<a) a=aaa
   enddo
   a=sqrt(a)  
   print*,"a",a
   if(a<7d-2)then;print*,"1***",a;a=0;endif
  ! if(a>0.006)a=a+a
   if(a>0.1)then;print*,"2***",a;a=a+a;endif
 !  if(a>0.02)a=a+a
   !EMD=sqrt(a)+EMD
   EMD=a+EMD
   print*,"EMD",EMD
  enddo
 ! a=huge(a)
  do i=1,nd
   a=1d4
   if(dese(i,1)>4000)cycle
   do j=1,ndop
    if(points0(j,1)>4000)cycle
    aaa=(dese(i,1)-points0(j,1))**2+(dese(i,2)-points0(j,2))**2+(dese(i,3)-points0(j,3))**2
 !   if(aaa<0.005)aaa=0
 !   if(aaa>0.005)aaa=aaa+aaa
    if(aaa<a) a=aaa
   enddo
   a=sqrt(a)  
   print*,"a",a
   if(a<7d-2)then;print*,"1***",a;a=0;endif
  ! if(a>0.006)a=a+a
   if(a>0.1)then;print*,"2***",a;a=a+a+a+a;endif
 !  if(a>0.02)a=a+a
   !EMD=sqrt(a)+EMD
   EMD=a+EMD
   print*,"EMD",EMD
  enddo
  EMD=EMD/real(ndtt+nddd)
  !!!print*,"EMD",EMD
print*,"emd",emd
  sce=EMD
  
   
!    write(39,*) pasadas,EMD

 

end subroutine compare2
!**************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compare3(dese,sce)
real*8 :: sce,aaa,a,EMD
!integer :: ndop
real*8 , dimension(nd,3) :: dese
 character*10 pasadach
!integer*8 :: ndtt
sce=0
!a=huge(a)
EMD=0.0
pasadas=pasadas+1
!ndop=ndtt
!write(pasadach,'(I10)')pasadas
  
! open(38,file=trim(descensi)//"_final_"//trim(adjustl(pasadach)))
!   do i=1,nd
!    write(38,*) dese(i,1),dese(i,2),dese(i,3)
!   enddo
! close(38)

  
  do i=1,ndop
   a=1d4
   if(points0(i,1)>4000)cycle
   do j=1,nd!dd
    if(dese(j,1)>4000)cycle
    aaa=(points0(i,1)-dese(j,1))**2+(points0(i,2)-dese(j,2))**2+(points0(i,3)-dese(j,3))**2
    if(aaa.ne.aaa)aaa=0
    if(aaa<a) a=aaa
   enddo
   EMD=sqrt(a)+EMD
  enddo
 ! a=huge(a)
  do i=1,nd
   a=1d4
   if(dese(i,1)>4000)cycle
   do j=1,ndop
    if(points0(j,1)>4000)cycle
    aaa=(dese(i,1)-points0(j,1))**2+(dese(i,2)-points0(j,2))**2+(dese(i,3)-points0(j,3))**2
    if(aaa.ne.aaa)aaa=0
    if(aaa<a) a=aaa
   enddo
   EMD=sqrt(a)+EMD
  enddo
  EMD=EMD/real(ndtt+nddd)
  !!!print*,"EMD",EMD

  sce=EMD
  
   
!    write(39,*) pasadas,EMD

 

end subroutine compare3
!**************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine compare(dese,sce)
real*8 :: sce,aaa,a,EMD,bbb,ccc
real*8 :: u1x,u1y,u1z,p1x,p1y,p1z,uxpx,uxpy,uxpz,Npux,Npuy,Npuz,nnv1,nnp
real*8 :: tnnx,tnny,tnnz,NN,t,p0x,p0y,p0z,uxpx3,mini
real*8 :: vx,vz,vy,sx,sy,sz,dopt
real*8 , dimension(nd,3) :: dese
integer , dimension(3) :: desetr
integer:: ff
 character*10 pasadach
integer*8 :: desetr1,desetr2,desetr3,ndoptimum,ndother
real*8 , allocatable :: aaa2(:)
sce=0
EMD=0.0

!find 3 closest points**
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(nd))
aaa2=1d4
do i=1,ndop
 if(points0(i,1)>4000)cycle
   do j=1,nd
    if(dese(j,1)>4000)cycle    
    aaa2(j)=(((points0(i,1)-dese(j,1))**2)+((points0(i,2)-dese(j,2))**2)+&
        &((points0(i,3)-dese(j,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,nd  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo


!vectors of triangle; normal vector of triangle
  
  u1x=dese(desetr(2),1)-dese(desetr(1),1);u1y=dese(desetr(2),2)-dese(desetr(1),2);u1z=dese(desetr(2),3)-dese(desetr(1),3)!v1
  p1x=dese(desetr(3),1)-dese(desetr(1),1);p1y=dese(desetr(3),2)-dese(desetr(1),2);p1z=dese(desetr(3),3)-dese(desetr(1),3)!v2
!  q1x=dese(desetr(3),1)-dese(desetr(1),1);q1y=dese(desetr(3),2)-dese(desetr(1),2);q1z=dese(desetr(3),3)-dese(desetr(1),3)!v3
  ptx=points0(i,1);pty=points0(i,2);ptz=points0(i,3)	!point
 
!unit normal vector for triangle

  uxpx3=u1y*p1z-u1z*p1y; uxpy3=u1z*p1x-u1x*p1z; uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3**2+uxpy3**2+uxpz3**2)
!print*,"uxpx3",uxpx3,"uxpy3",uxpy3,"uxpz3",uxpz3
!print*,"NN",NN
  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN
!print*,"npux",npux,"npuy",npux,"npuz",npuz
  !3
  ! vx=ptx-u1x;vy=pty-u1y;vz=ptz-u1z
    vx=ptx-dese(desetr(2),1);vy=pty-dese(desetr(2),2);vz=ptz-dese(desetr(2),3)
!print*,"vx",vx,"vy",vy,"vz",vz
  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz
!print*,"dotp",dotp
  !5 point to plane
   !sx=vx-Npux*dotp
   !sy=vy-Npuy*dotp
   !sz=vz-Npuz*dotp
   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp    
! print*,"sx",sx,"sy",sy,"sz",sz   
  !6 distance

 !  if(a<7d-)then;print*,"1***",a;a=0;endif
 !  if(a>0.1)then;print*,"2***",a;a=a+a+a+a;endif

   EMD=sqrt(sx*sx+sy*sy+sz*sz)+EMD 

!print*,dd!;stop

enddo
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(ndop))
aaa2=1d4
do i=1,nd
  if(dese(i,1)>4000)cycle
  do j=1,ndop
   if(points0(j,1)>4000)cycle
    aaa2(j)=(((points0(j,1)-dese(i,1))**2)+((points0(j,2)-dese(i,2))**2)+&
        &((points0(j,3)-dese(i,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,ndop  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo


 
  !vectors of triangle; normal vector of plane
  
  u1x=points0(desetr(2),1)-points0(desetr(1),1); u1y=points0(desetr(2),2)-points0(desetr(1),2)!v1
  u1z=points0(desetr(2),3)-points0(desetr(1),3)
  p1x=points0(desetr(3),1)-points0(desetr(1),1); p1y=points0(desetr(3),2)-points0(desetr(1),2)!v2
  p1z=points0(desetr(3),3)-points0(desetr(1),3)
!  q1x=points0(desetr(3),1)-points0(desetr(1),1);q1y=points0(desetr(3),2)-points0(desetr(1),2)!v3
!  q1z=points0(desetr(3),3)-points0(desetr(1),3)
  ptx=dese(i,1);pty=dese(i,2);ptz=dese(i,3)	!point
!normal vector
  uxpx3=u1y*p1z-u1z*p1y;uxpy3=u1z*p1x-u1x*p1z;uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3**2+uxpy3**2+uxpz3**2)

  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN

  !3
!   vx=ptx-u1x;vy=pty-u1y;vz=ptz-u1z
     vx=ptx-points0(desetr(2),1);vy=pty-points0(desetr(2),2);vz=ptz-points0(desetr(2),3)
  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz
  !5 point to plane
   !sx=vx-Npux*dotp
   !sy=vy-Npuy*dotp
   !sz=vz-Npuz*dotp
   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp   
  !6 distance
   EMD=sqrt(sx*sx+sy*sy+sz*sz)+EMD 
!   open(332,file="sxyz.dat")
!   write(332,*)sx,sy,sz
!   close(332)

!   open(333,file="Npuxyz.dat")
!   write(333,*)Npux,Npuy,Npuz
!   close(333)

!  open(334,file="ptxyz.dat")
!   write(334,*)ptx,pty,ptz
!   close(334)

!   open(335,file="u1xyz.dat")
!   write(335,*)u1x,u1y,u1z
!   close(335)
  
!   open(336,file="p1xyz.dat")
!   write(336,*)points0(desetr(2),1),points0(desetr(2),2),points0(desetr(2),3)
!   write(336,*)points0(desetr(1),1),points0(desetr(1),2),points0(desetr(1),3)
!   write(336,*)points0(desetr(3),1),points0(desetr(3),2),points0(desetr(3),3)
!   close(336)
!print*,EMD
!stop
enddo
!print*,dd,ndop,nd
!print*,dd/real(nddd+ndtt)!;stop
  !EMD=dd
  EMD=EMD/real(nddd+ndtt)
 ! EMD=EMD/real(ndoptimum+ndother)
!  print*,
 !   write(45,*)emd
  
  sce=EMD
! stop!ojo******** 
end subroutine compare

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!**************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine comparefinal(dese,sce)
real*8 :: sce,aaa,a,EMD,bbb,ccc,emd676
real*8 :: u1x,u1y,u1z,p1x,p1y,p1z,uxpx,uxpy,uxpz,Npux,Npuy,Npuz,nnv1,nnp
real*8 :: tnnx,tnny,tnnz,NN,t,p0x,p0y,p0z,uxpx3,mini
real*8 :: vx,vz,vy,sx,sy,sz,dopt
real*8 , dimension(nd,3) :: dese
integer, dimension(3) :: desetr
integer :: ff
 character*10 pasadach
integer*8 :: desetr1,desetr2,desetr3,ndoptimum,ndother
real*8 , allocatable :: aaa2(:)
sce=0
EMD=0.0

!find 3 closest points**
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(nd))
aaa2=1d4
do i=1,ndop
 if(points0(i,1)>4000)cycle
   do j=1,nd
    if(dese(j,1)>4000)cycle    
    aaa2(j)=(((points0(i,1)-dese(j,1))**2)+((points0(i,2)-dese(j,2))**2)+&
        &((points0(i,3)-dese(j,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,nd  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo


!vectors of triangle; normal vector of triangle
  
  u1x=dese(desetr(2),1)-dese(desetr(1),1);u1y=dese(desetr(2),2)-dese(desetr(1),2);u1z=dese(desetr(2),3)-dese(desetr(1),3)!v1
  p1x=dese(desetr(3),1)-dese(desetr(1),1);p1y=dese(desetr(3),2)-dese(desetr(1),2);p1z=dese(desetr(3),3)-dese(desetr(1),3)!v2
  ptx=points0(i,1);pty=points0(i,2);ptz=points0(i,3)	!point
 
!unit normal vector for triangle

  uxpx3=u1y*p1z-u1z*p1y; uxpy3=u1z*p1x-u1x*p1z; uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3**2+uxpy3**2+uxpz3**2)
  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN

  !3

    vx=ptx-dese(desetr(2),1);vy=pty-dese(desetr(2),2);vz=ptz-dese(desetr(2),3)

  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz

   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp    

  !6 distance
  emd676=sqrt(sx*sx+sy*sy+sz*sz)
 ! if(emd676>0.1)then
 !    emd676=emd676+emd676+emd676+emd676
 ! endif

   EMD=emd676*emd676+EMD 

enddo
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(ndop))
aaa2=1d4
do i=1,nd
  if(dese(i,1)>4000)cycle
  do j=1,ndop
   if(points0(j,1)>4000)cycle
    aaa2(j)=(((points0(j,1)-dese(i,1))**2)+((points0(j,2)-dese(i,2))**2)+&
        &((points0(j,3)-dese(i,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,ndop  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo


 
  !vectors of triangle; normal vector of plane
  
  u1x=points0(desetr(2),1)-points0(desetr(1),1); u1y=points0(desetr(2),2)-points0(desetr(1),2)!v1
  u1z=points0(desetr(2),3)-points0(desetr(1),3)
  p1x=points0(desetr(3),1)-points0(desetr(1),1); p1y=points0(desetr(3),2)-points0(desetr(1),2)!v2
  p1z=points0(desetr(3),3)-points0(desetr(1),3)

  ptx=dese(i,1);pty=dese(i,2);ptz=dese(i,3)	!point
!normal vector
  uxpx3=u1y*p1z-u1z*p1y;uxpy3=u1z*p1x-u1x*p1z;uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3**2+uxpy3**2+uxpz3**2)

  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN

  !3

     vx=ptx-points0(desetr(2),1);vy=pty-points0(desetr(2),2);vz=ptz-points0(desetr(2),3)
  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz
  !5 point to plane

   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp   
  !6 distance
!  emd676=sqrt(sx*sx+sy*sy+sz*sz)
!  if(emd676>0.1)then
!     emd676=emd676+emd676+emd676+emd676
!  endif

   EMD=emd676*emd676+EMD


  ! EMD=sqrt(sx*sx+sy*sy+sz*sz)+EMD 

enddo

  EMD=EMD/real(nddd+ndtt)

  sce=EMD

end subroutine comparefinal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine compare44(dese,sce)
real*8 :: sce,aaa,a,EMD,bbb,ccc
real*8 :: u1x,u1y,u1z,p1x,p1y,p1z,uxpx,uxpy,uxpz,Npux,Npuy,Npuz,nnv1,nnp
real*8 :: tnnx,tnny,tnnz,NN,t,p0x,p0y,p0z,uxpx3,mini
real*8 :: vx,vz,vy,sx,sy,sz,dopt,emd676
real*8 , dimension(nd,3) :: dese
integer, dimension(3) :: desetr
integer :: ff
 character*10 pasadach
integer*8 :: desetr1,desetr2,desetr3,ndoptimum,ndother
real*8 , allocatable :: aaa2(:)

open(101,file="desetr.dat")
open(102,file="point0.dat")


sce=0
EMD=0.0

!find 3 closest points**
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(nd))
aaa2=1d4
do i=1,ndop
 if(points0(i,1)>4000)cycle
   do j=1,nd
    if(dese(j,1)>4000)cycle    
    aaa2(j)=(((points0(i,1)-dese(j,1))**2)+((points0(i,2)-dese(j,2))**2)+&
        &((points0(i,3)-dese(j,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,nd  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo



!vectors of triangle; normal vector of triangle
  
  u1x=dese(desetr(2),1)-dese(desetr(1),1);u1y=dese(desetr(2),2)-dese(desetr(1),2);u1z=dese(desetr(2),3)-dese(desetr(1),3)!v1
  p1x=dese(desetr(3),1)-dese(desetr(1),1);p1y=dese(desetr(3),2)-dese(desetr(1),2);p1z=dese(desetr(3),3)-dese(desetr(1),3)!v2

  ptx=points0(i,1);pty=points0(i,2);ptz=points0(i,3)	!point
 
!unit normal vector for triangle

  uxpx3=u1y*p1z-u1z*p1y; uxpy3=u1z*p1x-u1x*p1z; uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3*uxpx3+uxpy3*uxpy3+uxpz3*uxpz3)

  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN
  !3
  vx=ptx-dese(desetr(2),1);vy=pty-dese(desetr(2),2);vz=ptz-dese(desetr(2),3)
  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz
  !5 point to plane
   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp     
  !6 distance
   emd676=sqrt(sx*sx+sy*sy+sz*sz)
   if(emd676>0.1)then
     print*,emd676,"node",i
     write(102,*)ptx,pty,ptz
   endif
    !if(emd676>0.1)then
!    if(i==1)then
!     print*,emd676,"node",i
!     print*,"desetr",desetr
!     print*,sx,sy,sz
  !   write(102,*)ptx,pty,ptz!
  !
  !   write(101,*)dese(desetr(1),1),dese(desetr(1),2),dese(desetr(1),3)
  !   write(101,*)dese(desetr(2),1),dese(desetr(2),2),dese(desetr(2),3)
  !   write(101,*)dese(desetr(3),1),dese(desetr(3),2),dese(desetr(3),3)
  !   print*,"*******************************1"
!    endif
 !  if(a<7d-)then;print*,"1***",a;a=0;endif
 !  if(a>0.1)then;print*,"2***",a;a=a+a+a+a;endif

   EMD=sqrt(sx*sx+sy*sy+sz*sz)+EMD 
!print*,"emd",emd
!print*,dd!;stop

enddo
!print*,dd!/real(ndop+nd)
!stop!OJO****************************///
if(allocated(aaa2))deallocate(aaa2)
allocate(aaa2(ndop))
aaa2=1d4
do i=1,nd
  if(dese(i,1)>4000)cycle
  do j=1,ndop
   if(points0(j,1)>4000)cycle
    aaa2(j)=(((points0(j,1)-dese(i,1))**2)+((points0(j,2)-dese(i,2))**2)+&
        &((points0(j,3)-dese(i,3))**2))
  enddo
 desetr=0
 do ff=1,3
  mini=1d6
  do j=1,nd  
   if(desetr(1)==j.or.desetr(2)==j)cycle
   if(mini>aaa2(j))then;mini=aaa2(j);desetr(ff)=j;endif
  enddo 
 enddo

!print*,"desetr",desetr
!print*,"mini",mini
!stop
 
  !vectors of triangle; normal vector of plane
  
  u1x=points0(desetr(2),1)-points0(desetr(1),1); u1y=points0(desetr(2),2)-points0(desetr(1),2)!v1
  u1z=points0(desetr(2),3)-points0(desetr(1),3)
  p1x=points0(desetr(3),1)-points0(desetr(1),1); p1y=points0(desetr(3),2)-points0(desetr(1),2)!v2
  p1z=points0(desetr(3),3)-points0(desetr(1),3)
!  q1x=points0(desetr(3),1)-points0(desetr(1),1);q1y=points0(desetr(3),2)-points0(desetr(1),2)!v3
!  q1z=points0(desetr(3),3)-points0(desetr(1),3)
  ptx=dese(i,1);pty=dese(i,2);ptz=dese(i,3)	!point
!normal vector
  uxpx3=u1y*p1z-u1z*p1y;uxpy3=u1z*p1x-u1x*p1z;uxpz3=u1x*p1y-u1y*p1x

  NN=sqrt(uxpx3**2+uxpy3**2+uxpz3**2)

  Npux=uxpx3/NN					
  Npuy=uxpy3/NN
  Npuz=uxpz3/NN

  !3
   vx=ptx-points0(desetr(2),1);vy=pty-points0(desetr(2),2);vz=ptz-points0(desetr(2),3)
   
  !4- dot product
   dotp=vx*Npux+vy*Npuy+vz*Npuz
  !5 point to plane
   !sx=vx-Npux*dotp
   !sy=vy-Npuy*dotp
   !sz=vz-Npuz*dotp
   sx=-Npux*dotp
   sy=-Npuy*dotp
   sz=-Npuz*dotp   
  !6 distance
!    emd676=sqrt(sx*sx+sy*sy+sz*sz)
!    if(emd676>0.3)print*,emd676,"node",i
!print*,"*******************************2"
   EMD=sqrt(sx*sx+sy*sy+sz*sz)+EMD 
!   open(332,file="sxyz.dat")
!   write(332,*)sx,sy,sz
!   close(332)

!   open(333,file="Npuxyz.dat")
!   write(333,*)Npux,Npuy,Npuz
!   close(333)

!  open(334,file="ptxyz.dat")
!   write(334,*)ptx,pty,ptz
!   close(334)

!   open(335,file="u1xyz.dat")
!   write(335,*)u1x,u1y,u1z
!   close(335)
  
!   open(336,file="p1xyz.dat")
!   write(336,*)points0(desetr(2),1),points0(desetr(2),2),points0(desetr(2),3)
!   write(336,*)points0(desetr(1),1),points0(desetr(1),2),points0(desetr(1),3)
!   write(336,*)points0(desetr(3),1),points0(desetr(3),2),points0(desetr(3),3)
!   close(336)
!print*,EMD
!stop
enddo
!print*,dd,ndop,nd
!print*,dd/real(nddd+ndtt)!;stop
print*,"nddd",nddd,"ndtt",ndtt
print*,"ndoptimum",ndoptimum,"ndother",ndother
  !EMD=dd
  EMD=EMD/real(ndoptimum+ndother)
  
 !   write(45,*)emd
  
  sce=EMD
! stop!ojo******** 
 close(101);close(102)
end subroutine compare44



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine compare(boxes1,boxcount0,sc)

!integer :: mmx,mmy,mmz,nnx,nny,nnz,i,ii,iii
!integer, allocatable :: boxes1(:,:,:)
!integer, dimension(boxn) :: boxcount0
!real*8 :: sc,boxcountr

!boxcount=0; sc=0; 
!mmx=mx;mmy=my;mmz=mz;nnx=nx;nny=ny;nnz=nz

! this part is meant to make the code faster: Since in the majority of cases, the morphologies will fit into the same box number, 
! we pass through all the boxes1 separately only if the box numbers differ.

!if((mx.lt.step-1).or.(my.lt.step-1).or.(mz.lt.step-1).or.(nx.gt.0).or.(ny.gt.0).or.(nz.gt.0))then  !>>> Is 24-6-14
! there is some overlap between the boxes
! if(mx.le.0) mmx=0
! if(my.le.0) mmy=0
! if(mz.le.0) mmz=0
! if(nx.ge.step-1) nnx=step-1
! if(ny.ge.step-1) nny=step-1
! if(nz.ge.step-1) nnz=step-1

! do i=mmx,nnx
!   do ii=mmy,nny
!     do iii=mmz,nnz
!       if(boxes1(i,ii,iii).eq.1)then
!         boxcount=boxcount+1
!         if(boxes0(i,ii,iii,step).eq.1) sc=sc+1 ! matched boxes
!       endif
!     end do
!   end do
! end do

!end if
! if the box numbers differ, count again all filled boxes in the descendant (the ancestor's ones are already counted)
!if((mx.ne.0).or.(my.ne.0).or.(mz.ne.0).or.(nx.ne.step-1).or.(ny.ne.step-1).or.(nz.ne.step-1))then

! boxcount=0
! do i=mx,nx
!   do ii=my,ny
!     do iii=mz,nz
!       if(boxes1(i,ii,iii).eq.1) boxcount=boxcount+1  
!     end do
!   end do
! end do
!end if

!boxcount=boxcount+boxcount0(step) !>>> RZ 29-06
!boxcountr=boxcount
!sc=1.0d0-sc*2d0/boxcountr ! >>> RZ 29-06 we need to duplicate sc in order to account for both !!boxcounts

!end subroutine compare
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ANCESTOR SUBROUTINES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ancestor_define_boxes(ancestori,aao,boxnum2,sizedifi)

integer :: steps,nodj,nd0,boxnum2
real*8 ::  x,y,z,a,aao,sizedifi
real*8,dimension(3) :: medi
character*140 ancestori,ancestorii

call system("cp "//trim(ancestori)//" test00 ")

nd0=0
open(2,file=ancestori)
do j=1,1000000
  read(2,*,END=109,ERR=109) x, y, z, a, a, a
  nd0=nd0+1
end do
close(2)

109 print*, "READIN DONE", ancestori, " POINTS: ", nd0
if(allocated(points0)) deallocate(points0)
allocate(points0(nd0,3))

open(2,file="test00")
do j=1,nd0
  read(2,*,END=111,ERR=111) x, y, z, a, a, a
  points0(j,1)=x
  points0(j,2)=y
  points0(j,3)=z
end do
111 close(2)
nd=nd0
ndtt=nd
steps=1

!if(allocated(boxes0)) deallocate(boxes0)
!if(allocated(boxcount0)) deallocate(boxcount0)
!if(allocated(c0x)) deallocate(c0x)
!if(allocated(c0y)) deallocate(c0y)
!if(allocated(c0z)) deallocate(c0z)
!allocate(c0x(1),c0y(1),c0z(1))
!allocate(boxes0(0:boxnum2-1,0:boxnum2-1,0:boxnum2-1,boxnum2),boxcount0(boxnum2),c0x(boxnum2),c0y(boxnum2),c0z(boxnum2))
! the reason why I have deviated from the suggestion to center the boxes around 0,0,0 is that this raises problems for odd-numbered versus even-numberd grid sizes. Instead, I center around 0.5,0.5,0.5 times grid size.

 minxc=huge(minxc);minyc=huge(minyc);minzc=huge(minzc);maxxc=-1*huge(maxxc);maxyc=-1*huge(maxyc);maxzc=-1*huge(maxzc)
 maxvo=0.0d0
 c0x=0.0d0; c0y=0.0d0; c0z=0.0d0

 do j=1,nd ! save the values of the minimal nodes of x,y,z
   if(points0(j,1).lt.minxc) minxc=points0(j,1)
   if(points0(j,2).lt.minyc) minyc=points0(j,2)
   if(points0(j,3).lt.minzc) minzc=points0(j,3)
   if(points0(j,1).gt.maxxc) maxxc=points0(j,1)
   if(points0(j,2).gt.maxyc) maxyc=points0(j,2)
   if(points0(j,3).gt.maxzc) maxzc=points0(j,3)
 end do 

 maxvo=maxxc-minxc ! now search for the maximal extension in whatever direction; it is important for the absolute size
 if(maxvo.lt.maxyc-minyc) maxvo=maxyc-minyc
 if(maxvo.lt.maxzc-minzc) maxvo=maxzc-minzc
 boxes0=0; boxcount0=0
! aao=0d0; pp=0d0 pfh
 do nodj=1,nd
   call ancestor_fill_boxes(nodj); c0y=c0y+(points0(nodj,2))-miny; c0z=c0z+(points0(nodj,3))-minz  ! >>> RZ 29-06 centroid 
 end do
   !do j=1,nd
   !if(nod==j)cycle
  !   p=(points0(nod,1)-points0(j,1))**2+(points0(nod,2)-points0(j,2))**2+(points0(nod,3)-points0(j,3))**2
 !    if(p>pp) pp=p
 !  end do
! end do
! pp=sqrt(pp)
! sizedifi=pp
! print*, "PP", pp

 c0x=c0x/nd; c0y=c0y/nd; c0z=c0z/nd ! >>> RZ 29-06 centroid 

 !steps=boxnum2  
   !boxcount0(steps)=sum(boxes0(:,:,:,steps)) ! the total number of occupied boxes
   !c0x(steps)=c0x(steps)/maxvo*steps; c0y(steps)=c0y(steps)/maxvo*steps; c0z(steps)=c0z(steps)/maxvo*steps ! >> RZ 29-06 centroid
 

end subroutine ancestor_define_boxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ancestor_fill_boxes(nodj)!,ooa,nd00,boxnum3)!!

integer :: nodj,steps,bx,by,bz,nd00,boxnum3
real*8 :: ax,ay,az,ooa
! shift all nodes so that they range between 0 and maxvo and that their centre is in 0.5,0.5,0.5 times maxvo

   ax=points0(nodj,1)-minxc
   ay=points0(nodj,2)-minyc
   az=points0(nodj,3)-minzc
   ax=ax/maxvo ! --> all values range between 0 and 1
   ay=ay/maxvo
   az=az/maxvo

 !  do j=1,nd00
!     b=sqrt((ax-(points0(j,1)-minx)/maxvo)**2+(ay-(points0(j,2)-miny)/maxvo)**2+(az-(points0(j,3)-minz)/maxvo)**2)
!     if(b>ooa) ooa=b
!   end do

!if(nod==nd00) print*, "HEUREKA", ooa
open(17,file="AA",position="append")
   !do steps=boxs,boxn
   steps=1
     bx=int(steps*ax); by=int(steps*ay); bz=int(steps*az) ! fill into boxes
     if((points0(nodj,1)-minxc)==maxvo) bx=bx-1 ! this is for the extreme node
     if((points0(nodj,2)-minyc)==maxvo) by=by-1
     if((points0(nodj,3)-minzc)==maxvo) bz=bz-1
     boxes0(bx,by,bz,steps)=1
     write(17,*) bx, by, bz
   !end do
close(17)
end subroutine ancestor_fill_boxes



end module
