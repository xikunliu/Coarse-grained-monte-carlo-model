
subroutine CalcEnergy(Pos, AA, M, L, rc, BL, SL, sigma, PEnergy, Dim, NAtom)
    implicit none
    integer, intent(in) :: M, Dim, NAtom
    real(8), intent(in), dimension(0:NAtom-1, 0:Dim-1) :: Pos
    real(8), intent(in), dimension(0:NAtom-1) :: AA
    real(8), intent(in) :: L, rc
    real(8), intent(in) :: BL,SL,sigma
    real(8), intent(out) :: PEnergy
    real(8), parameter :: k = 236., r0 = 3.8
    !======== YOUR CODE HERE ========
    real(8), dimension(Dim) :: rij, Posi
    real(8) :: d, d2, id2, sid2, sid6, sid12
    real(8) ::  sigma2
    real(8) :: rc2, irc2, sirc2, sirc6, sirc12, Shift1,Shift2
    integer :: i, j
    PEnergy = 0.
    rc2 = rc * rc
    sigma2 = sigma * sigma      
    irc2 = 1. /rc2
    sirc2 = sigma2 * irc2
    sirc6 = sirc2*sirc2*sirc2
    sirc12=sirc6*sirc6
    Shift1 = -4. * 5.* (sirc12 - sirc6)
    Shift2 = -4.* (sirc12 - sirc6)
    do i = 0, NAtom - 1
        !store Pos(i,:) in a temporary array for faster access in j loop
        Posi = Pos(i,:)
        do j = i + 1, NAtom - 1
            rij = Pos(j,:) - Posi
            rij = rij - L * dnint(rij / L)
            !compute only the squared distance and compare to squared cut
            d2 = sum(rij * rij)
            if (d2 > rc2) then
                cycle
            end if
            if (j==i+1 .and. mod(j,M)>0)  then !covalent
                d=sqrt(d2)
                PEnergy = PEnergy + k/2*(d-r0)**2
            ! coulomb interation with same charges
            else if ((AA (i) == 1 .and. AA (j) == 1) .or. (AA (i) == 2 .and. AA (j) == 2)) then
                d=sqrt(d2)
                PEnergy = PEnergy + BL/d * exp (-d/SL)
            ! coulomb interation with opposite charges
            else if ((AA (i) == 1 .and. AA (j) == 2) .or. (AA (i) == 2 .and. AA (j) == 1)) then
                d=sqrt(d2)
                PEnergy = PEnergy - BL/d * exp (-d/SL)
            ! interaction between two hydrophobic residues
            else if (AA(i) == 3 .and. AA(j) == 3) then
                sigma2 = sigma * sigma
                id2 = 1. / d2            !inverse squared distance 
                sid2=sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma       
                PEnergy = PEnergy + 4.* 5. * (sid12 - sid6) + Shift1
            ! interaction between two hydrophilic residues or different residues
            else
                sigma2 = sigma * sigma
                id2 = 1. / d2            !inverse squared distance 
                sid2=sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma       
                PEnergy = PEnergy + 4. * (sid12 - sid6) + Shift2
            endif
        enddo
    enddo
end subroutine

subroutine CalcEnergySpecified(Pos, AA, M, L, rc, BL, SL, sigma, i, PEnergy, Dim, NAtom)
    implicit none
    integer, intent(in) :: M, i, Dim, NAtom
    real(8), intent(in), dimension(0:NAtom-1, 0:Dim-1) :: Pos
    real(8), intent(in), dimension(0:NAtom-1) :: AA
    real(8), intent(in) :: L, rc, BL, SL, sigma
    real(8), intent(out) :: PEnergy
    real(8), parameter :: k = 236., r0 = 3.8
    !======== YOUR CODE HERE ========
    real(8), dimension(Dim) :: rij, Posi
    real(8) :: d, d2, id2, sid2, sid6, sid12
    real(8) :: sigma2
    real(8) :: rc2, irc2, sirc2, sirc6, sirc12, Shift1,Shift2
    integer :: j
    PEnergy = 0.
    rc2 = rc * rc
    sigma2 = sigma * sigma      
    irc2 = 1. /rc2
    sirc2 = sigma2 * irc2
    sirc6 = sirc2*sirc2*sirc2
    sirc12=sirc6*sirc6
    Shift1 = -4. * 5.* (sirc12 - sirc6)
    Shift2 = -4.* (sirc12 - sirc6)

        !store Pos(i,:) in a temporary array for faster access in j loop
        Posi = Pos(i,:)
        do j = 0, NAtom - 1
            if (j == i) then
                cycle
            end if
            rij = Pos(j,:) - Posi
            rij = rij - L * dnint(rij / L)
            !compute only the squared distance and compare to squared cut
            d2 = sum(rij * rij)
            if (d2 > rc2) then
                cycle
            end if
            if ((j==i+1 .and. mod(j,M)>0) .or. (j==i-1 .and. mod(i,M)>0))  then !covalent
                d=sqrt(d2)
                PEnergy = PEnergy + k/2*(d-r0)**2
            ! coulomb interation with same charges
            else if ((AA (i) == 1 .and. AA (j) == 1) .or. (AA (i) == 2 .and. AA (j) == 2)) then
                d=sqrt(d2)
                PEnergy = PEnergy + BL/d * exp (-d/SL)
            ! coulomb interation with opposite charges
            else if ((AA (i) == 1 .and. AA (j) == 2) .or. (AA (i) == 2 .and. AA (j) == 1)) then
                d=sqrt(d2)
                PEnergy = PEnergy - BL/d * exp (-d/SL)
            ! interaction between two hydrophobic residues
            else if (AA(i) == 3 .and. AA(j) == 3) then
                id2 = 1. / d2            !inverse squared distance 
                sid2=sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma       
                PEnergy = PEnergy + 4.* 5. * (sid12 - sid6) + Shift1
            ! interaction between two hydrophilic residues or different residues
            else   
                id2 = 1. / d2            !inverse squared distance 
                sid2=sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma                
                PEnergy = PEnergy + 4. * (sid12 - sid6) + Shift2
            endif
        enddo

end subroutine

subroutine CalcEnergyForces(Pos, AA, M, L, rc, BL, SL, sigma, PEnergy, Forces, Dim, NAtom)
    implicit none
    integer, intent(in) :: M, Dim, NAtom
    real(8), intent(in), dimension(0:NAtom-1, 0:Dim-1) :: Pos
    real(8), intent(in), dimension(0:NAtom-1) :: AA
    real(8), intent(in) :: L, rc
    real(8), intent(in) :: BL,SL,sigma
    real(8), intent(out) :: PEnergy
    real(8), intent(inout), dimension(0:NAtom-1, 0:Dim-1) :: Forces
!f2py intent(in, out) :: Forces
    real(8), parameter :: k = 236., r0 = 3.8
    !======== YOUR CODE HERE ========
    real(8), dimension(Dim) :: rij, Fij, Posi
    real(8) :: d, id, d2, id2, sid2, sid6, sid12, Shift1, Shift2,rc2, irc2, sirc2, sirc6, sirc12
    real(8) ::  sigma2
    integer :: i, j
    PEnergy = 0.
    Forces = 0.
    rc2 = rc * rc
    sigma2 = sigma * sigma      
    irc2 = 1. /rc2
    sirc2 = sigma2 * irc2
    sirc6 = sirc2*sirc2*sirc2
    sirc12=sirc6*sirc6
    Shift1 = -4. * 5.* (sirc12 - sirc6)
    Shift2 = -4.* (sirc12 - sirc6)
    do i = 0, NAtom - 1
        !store Pos(i,:) in a temporary array for faster access in j loop
        Posi = Pos(i,:)
        do j = i + 1, NAtom - 1
            rij = Pos(j,:) - Posi
            rij = rij - L * dnint(rij / L)
            !compute only the squared distance and compare to squared cut
            d2 = sum(rij * rij)
            if (d2 > rc2) then
                cycle
            end if
        !inverse twelvth distance
            if (j==i+1 .and. mod(j,M)>0)  then
                d=sqrt(d2)
                id=1. /d
                PEnergy = PEnergy + k/2*(d-r0)**2
                Fij = rij*k*id*(d-r0)
                Forces(i,:) = Forces(i,:) + Fij
                Forces(j,:) = Forces(j,:) - Fij
            else if ((AA (i) == 1 .and. AA (j) == 1) .or. (AA (i) == 2 .and. AA (j) == 2)) then
                d=sqrt(d2)
                id = 1. /d
                id2=id*id
                PEnergy = PEnergy + BL*id * exp (-d/SL)
                Fij = -rij*id*(BL*id2+BL*id*(1/SL))* exp(-d/SL)
                Forces(i,:) = Forces(i,:) + Fij
                Forces(j,:) = Forces(j,:) - Fij
            else if ((AA (i) == 1 .and. AA (j) == 2) .or. (AA (i) == 2 .and. AA (j) == 1)) then
                d=sqrt(d2)
                id = 1/d
                id2=id*id
                PEnergy = PEnergy + BL*id * exp (-d/SL)
                Fij = rij*id*(BL*id2+BL*id*(1/SL))* exp(-d/SL)
                Forces(i,:) = Forces(i,:) + Fij
                Forces(j,:) = Forces(j,:) - Fij
            ! interaction between two hydrophobic residues
            else if (AA(i) == 3 .and. AA(j) == 3) then
                id2 = 1. / d2            !inverse squared distance 
                sid2=sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma       
                PEnergy = PEnergy + 4.* 5. * (sid12 - sid6) + Shift1
                Fij = rij * ((-240. * sid12 + 120. * sid6) * sid2)                              
            else
                id2 = 1. / d2            !inverse squared distance 
                sid2 = sigma2*id2          ! inverse squared distance with sigma
                sid6 = sid2 * sid2 * sid2    !inverse sixth distance with sigma
                sid12 = sid6 * sid6         !inverse twelvth distance with sigma       
                PEnergy = PEnergy +  4. * (sid12 - sid6) + Shift2
                Fij = rij * ((-48. * sid12 + 24. * sid6) * sid2) 
            end if
        enddo
    enddo
end subroutine



subroutine CalcContacts(Pos, L, Cut, Contacts, Avg, Dim, NAtom)
    implicit none
    integer, intent(in) :: Dim, NAtom
    real(8), intent(in), dimension(0:NAtom-1, 0:Dim-1) :: Pos
    real(8), intent(in) :: Cut,L
    integer, intent(out), dimension(0:NAtom-1) :: Contacts
    real(8), intent(out) :: Avg
    real(8), dimension(Dim) :: rij, Posi
    real(8) :: d2
    real(8) :: Cut2
    integer :: i, j
    Contacts = 0
    Cut2 = Cut * Cut
    do i = 0, NAtom - 1
        !store Pos(i,:) in a temporary array for faster access in j loop
        Posi = Pos(i,:)
        do j = i + 1, NAtom - 1
            rij = Pos(j,:) - Posi
            rij = rij - L * dnint(rij / L)
            !compute only the squared distance and compare to squared cut
            d2 = sum(rij * rij)
            if (d2 < Cut2) then
                Contacts(i)= Contacts(i)+1
                Contacts(j)=Contacts(j)+1
            end if
        enddo
    enddo
    Avg = sum(Contacts)/NAtom
end subroutine

subroutine CalcContactsSpecified(Nchain, Contacts, NResidue, Avgs, NAtom)
    implicit none
    integer, intent(in) :: NAtom, Nchain
    integer, intent(in) :: NResidue
    integer, intent(in), dimension(0:NAtom-1) :: Contacts
    real(8), intent(out) :: Avgs
    integer :: k
    integer :: sums
    sums = 0
    do k = 0, NAtom - 1
        if (mod(k,8) == NResidue) then
            sums=sums+Contacts(k)
        end if
    enddo
    Avgs = sums/Nchain 
end subroutine

subroutine CalcContactsAA(AA, Contacts, NAA, AvgAA, NAtom)
    implicit none
    integer, intent(in) :: NAtom, NAA
    real(8), intent(in), dimension(0:NAtom-1) :: AA
    integer, intent(in), dimension(0:NAtom-1) :: Contacts
    real(8), intent(out) :: AvgAA
    integer :: k
    integer :: sumAA
    sumAA = 0 ! sum of contact number
    do k = 0, NAtom - 1
        if (AA(k) == NAA) then
            sumAA=sumAA+Contacts(k)
        end if
    enddo
    AvgAA = sumAA
end subroutine