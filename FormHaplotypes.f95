!! Written By David Howard
!! 09-10-2015
!! Last edited 01-11-2015

!! Takes phased data and builds haplotype library and outputs those above a specified MAF
!! Requires genetic map and .haps files (for format see SHAPEIT2) and uses a .fam file to find number of individuals
!! On command line, three values required: window size (cM), minimum number of markers in a window and minor haplotype frequency"

PROGRAM FormHaplotypes
IMPLICIT NONE
CHARACTER :: cn
CHARACTER(LEN=2) :: chromosome
CHARACTER(LEN=8), DIMENSION(:), ALLOCATABLE :: arg
CHARACTER(LEN=500), DIMENSION(:), ALLOCATABLE :: workhap, haplib
CHARACTER(LEN=4), DIMENSION(:,:), ALLOCATABLE :: output
CHARACTER(LEN=10), DIMENSION(:,:), ALLOCATABLE :: gmap, haplo, geno, ped
INTEGER :: ichr, ip, iq, ir=0, iu, iv, iw, ix, iy, iz, np=0, nq=0, nr=0, ierr, posstart, posend, nextwin, hapfound, obt, minwin
INTEGER, DIMENSION(:), ALLOCATABLE :: hapcount
INTEGER, PARAMETER :: dp = kind(1.d0)
REAL(kind=dp) :: wstart, maf, wsize
REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: cMpos

allocate(arg(3))
IF (IARGC() < 3) THEN
  print*, "3 command line arguments required: window size (cM), minimum number of markers in a window and minor haplotype frequency"
  STOP
ENDIF
CALL get_command_argument(1, arg(1))
read(arg(1),*) wsize
CALL get_command_argument(2, arg(2))
read(arg(2),*) minwin
CALL get_command_argument(3, arg(3))
read(arg(3),*) maf

print*, "Window size (cM) :", wsize
print*, "Minimum markers in a window :", minwin
print*, "Allowable minor haplotype frequency :", maf

open (unit=40,file="mapping.txt", STATUS="REPLACE")  !! File providing Mb to cM position
write(40,*), "Chromosome Mb cM"

open (unit=50,file="windows_used_"//trim(arg(1))//"_"//trim(arg(2))//"_"//trim(arg(3))//".txt", STATUS="REPLACE")
write(50,*), "Chromosome, Start_SNP End_SNP Start_cM End_cM #SNPsInWindow TotalHaplotypes QualifiedHaplotypes"

open (unit=60,file="haplotypes_"//trim(arg(1))//"_"//trim(arg(2))//"_"//trim(arg(3))//".map", STATUS="REPLACE")

open (unit=70,file="haplotypes_"//trim(arg(1))//"_"//trim(arg(2))//"_"//trim(arg(3))//".ped", STATUS="REPLACE")

Open (unit=10, file="QCdGS20K_chr1.fam", action='read', status="old")  ! Load input
LoadPedFile: DO
  READ(10,'(a)', IOSTAT=ierr) cn
  IF (ierr /= 0) EXIT
  nr = nr + 1		! Find total number of rows and set as array length nr
END DO LoadPedFile
Close (unit=10)

allocate(ped(6,nr))
Open (unit=11, file="QCdGS20K_chr1.fam", status="old")  ! Load input
Read(11,*) ped
Close (unit=11)

print*, nr, "individuals"
allocate(output(1000000,nr))  !! Set large default number of snps
output = ""

EachChr: do ichr = 1,22

np = 0
nq = 0

write (chromosome,'(i2)') ichr
chromosome = adjustl(chromosome)

Open (unit=10, file="genetic_map_chr"//trim(chromosome)//"_combined_b37.txt", action='read', status="old")  ! Load input
LoadMapFile: DO
  READ(10,'(a)', IOSTAT=ierr) cn
  IF (ierr /= 0) EXIT
  np = np + 1		! Find total number of rows and set as array length np
END DO LoadMapFile
Close (unit=10)
np = np - 1 !! Remove header row

Open (unit=10, file="chr"//trim(chromosome)//".phased.haps", action='read', status="old")  ! Load input
LoadHapsFile: DO
  READ(10,'(a)', IOSTAT=ierr) cn
  IF (ierr /= 0) EXIT
  nq = nq + 1		! Find total number of rows and set as array length nq
END DO LoadHapsFile
Close (unit=10)

print*, ""
print*, "Chromosome", ichr
print*, np, "markers in genetic map"
print*, nq, "markers in phased data"

allocate(workhap(nr*2+5))
allocate(haplib(nr))
allocate(hapcount(nr))

allocate(gmap(3,np))
Open (unit=11, file="genetic_map_chr"//trim(chromosome)//"_combined_b37.txt", status="old")  ! Load input
Read(11,*)
Read(11,*) gmap
Close (unit=11)

allocate(cMpos(nq))
cMpos = 0

allocate(haplo((nr*2)+5,nq))
Open (unit=11, file="chr"//trim(chromosome)//".phased.haps", action="read", status="old")  ! Load input
Read(11,*) haplo
Close (unit=11)

FindcM: do ix = 1,nq                     !! For each phase position
  ip = 0
  FromGMap: do iy = 1,np                 !! Look within genetic map
    if (haplo(3,ix) == gmap(1,iy)) then  !! If phase position in Mb is found in genetic map
	  read(gmap(3,iy),*) cMpos(ix)       !! Store cM position as real number within cMpos
	  write(40,*), ichr, haplo(3,ix), gmap(3,iy)
	  ip = 1
	  exit FromGMap
	endif
  enddo FromGMap
  if (ip == 0) then
    write(40,*), ichr, haplo(3,ix), "0"
  endif
enddo FindcM

wstart = cMpos(1)  !! starting cM position of a window
posstart = 1	   !! starting row position of a window
iq = 0   !! Number of haplotypes in each window
ix = 0   !! while loop below
Markers: do while (ix < nq)       !! For each marker position
  ix = ix + 1
  if ((cMpos(ix) /= 0) .AND. (cMpos(ix) <= wstart+(wsize/2))) then  !! Find next position of offset window
    nextwin = ix + 1
  endif
  if ((cMpos(ix) >= wstart+wsize) .OR. (ix == nq)) then  !! if phase is above the defined window length or at end of chromosome
	posend = ix-1                      !! posend takes value of previous row
	if (posend-posstart+1 >= minwin) then   !! if there are more than minwin markers in window
	  workhap = ""
	  if (posend-posstart+1 > 500) then !! if there are more than 500 phases in window need to amend code
	    print*, "Need to extend number of characters assigned to workhap and haplib to ", posend-posstart+1
	    STOP
	  endif
	  haplib = ""   !! Array of haplotypes observed
	  hapcount = 0  !! Array of number of times each haplotype observed
	  hapfound = 0  !! Number of unique haplotypes
	  ForEachInd: do iy = 6,(nr*2+5)  !! For each individual's phase
		obt = 0
		ObtainHap: do iz = posstart,posend   !! Use windows to define haplotype
		  obt = obt + 1
		  workhap(iy)(obt:obt) = haplo(iy,iz)    !! Temporarily store each haplotype in workhap
		enddo ObtainHap
!		print*, posend, obt, workhap(1:obt)(iy)
		if (hapfound > 0) then !! If not first haplotype
		  ip = 0 !! Haplotype observation
		  CheckLib: do iz = 1,hapfound
		    if (workhap(iy) == haplib(iz)) then  !! Haplotype already observed
			  hapcount(iz) = hapcount(iz) + 1
			  ip = 1
			endif
		  enddo CheckLib
		  if (ip == 0) then !! Haplotype not already observed
		    hapfound = hapfound + 1
			haplib(hapfound) = workhap(iy)
			hapcount(hapfound) = hapcount(hapfound) + 1
		  endif
		else !! First haplotype
	      hapfound = hapfound + 1
          haplib(hapfound) = workhap(iy)
		  hapcount(hapfound) = hapcount(hapfound) + 1
		endif
	  enddo ForEachInd

	  ForEachHap: do iz = 1,hapfound
	    if (hapcount(iz) > 2*nr*maf) then  !! Observed haplotypes with a maf above pre-defined value
		  ir = ir + 1
		    if (ir > 1000000) then
              print*, "Need to amend code for output matrix size"
            endif
		  iq = iq + 1
		  write(60,*), chromosome, ir, "0 ", haplo(3,posstart)
		  IdentifyInd: do iu = 6,(nr*2+5)  !! Prepare output for each individual
		    if (workhap(iu)(1:obt) == haplib(iz)(1:obt)) then  !! If individual carries haplotype of interest
			  if (MOD(iu,2) == 0) then !! if even and therefore first position
				output(ir,(iu/2)-2)(2:2) = "C"  !! Record as C
		      else !! odd and therefore in second position
			    output(ir,((iu-1)/2)-2)(4:4) = "C"  !! Record as C
			  endif
			else  !! If individual carries wild type haplotype
			  if (MOD(iu,2) == 0) then !! if even
			    output(ir,(iu/2)-2)(2:2) = "A"
		      else !! odd
			    output(ir,((iu-1)/2)-2)(4:4) = "A"
			  endif
			endif
		  enddo IdentifyInd

		end if
	  enddo ForEachHap

	write(50,*), ichr, posstart, posend, cMpos(posstart), cMpos(posend), posend-posstart+1, hapfound, iq
	endif

	iq = 0
	ix = nextwin
	posstart = nextwin
	wstart = wstart+(wsize/2)

  endif
enddo Markers

deallocate(workhap)
deallocate(haplib)
deallocate(hapcount)
deallocate(gmap)
deallocate(cMpos)
deallocate(haplo)

enddo EachChr

print*, ""
print*, ir, "total qualified haplotypes"

outputfile: do ix = 1,nr
  write(70,*), trim(ped(1,ix)), " ", trim(ped(2,ix)), " ", trim(ped(3,ix)), " ", trim(ped(4,ix)), " ", &
  trim(ped(5,ix)), " ", trim(ped(6,ix)), output(1:ir,ix)
enddo outputfile


END PROGRAM FormHaplotypes
