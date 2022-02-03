! Wrapper for initialization routine
!
! For documentation, see rrtmg_sw_init.f90
!
! Modification w/r/to rrtmg_sw_init:
! - change real(kind=rb) to real(8)
!   (f2py has trouble with parameterized types)
subroutine rrtmg_sw_ini_wrapper(cpdair)
  
  use rrtmg_sw_init, only: rrtmg_sw_ini
  implicit none
  real(8), intent(in) :: cpdair

  call rrtmg_sw_ini(cpdair)

end subroutine rrtmg_sw_ini_wrapper

! Wrapper for radiative transfer calculation
!
! For documentation, see rrtmg_sw_rad.f90
!
! Modification w/r/to rrtmg_sw_rad:
! - change icld and iaer from intent(inout) to intent(out)
!   to clarify that modifications will not be reflected
!   in value of input variable after subroutine returns
!   (avoid different behaviors for Python scalars and rank-0 arrays)
! - change real(kind=rb) to real(8)
!   and integer(kind=im) to integer(4)
!   (f2py has trouble with parameterized types)
! - convert intent(out) arrays to intent(inout)
!   (f2py has trouble with intent(inout) assumed shape arrays)
! - convert optional arguments to required arguments
!   (not sure how f2py handles optional arguments)
subroutine rrtmg_sw_wrapper(                        &
    ncol, nlay, icld, iaer,                         &
    play, plev, tlay, tlev, tsfc,                   &
    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,   &
    asdir, asdif, aldir, aldif,                     &
    coszen, adjes, dyofyr, scon, isolvar,           &
    inflgsw, iceflgsw, liqflgsw, cldfr,             &
    taucld, ssacld, asmcld, fsfcld,                 &
    cicewp, cliqwp, reice, reliq,                   &
    tauaer, ssaaer, asmaer, ecaer,                  &
    start_band, end_band,                           &
    swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc,  &
    wavenumber_range,                               &
    bndsolvar, indsolvar, solcycfrac )

    use rrtmg_sw_rad, only: rrtmg_sw
    implicit none
      
    integer(4), intent(in) :: ncol            
    integer(4), intent(in) :: nlay            
    integer(4), intent(in) :: icld         
    integer(4), intent(in) :: iaer         
    real(8), intent(in) :: play(:,:)          
    real(8), intent(in) :: plev(:,:)          
    real(8), intent(in) :: tlay(:,:)          
    real(8), intent(in) :: tlev(:,:)          
    real(8), intent(in) :: tsfc(:)            
    real(8), intent(in) :: h2ovmr(:,:)        
    real(8), intent(in) :: o3vmr(:,:)         
    real(8), intent(in) :: co2vmr(:,:)        
    real(8), intent(in) :: ch4vmr(:,:)        
    real(8), intent(in) :: n2ovmr(:,:)        
    real(8), intent(in) :: o2vmr(:,:)         
    real(8), intent(in) :: asdir(:)           
    real(8), intent(in) :: aldir(:)           
    real(8), intent(in) :: asdif(:)           
    real(8), intent(in) :: aldif(:)           
    integer(4), intent(in) :: dyofyr          
    real(8), intent(in) :: adjes              
    real(8), intent(in) :: coszen(:)          
    real(8), intent(in) :: scon               
    integer(4), intent(in) :: isolvar         
    real(8), intent(in) :: indsolvar(:) 
    real(8), intent(in) :: bndsolvar(:)  
    real(8), intent(in) :: solcycfrac    
    integer(4), intent(in) :: inflgsw        
    integer(4), intent(in) :: iceflgsw       
    integer(4), intent(in) :: liqflgsw       
    real(8), intent(in) :: cldfr(:,:)        
    real(8), intent(in) :: taucld(:,:,:)     
    real(8), intent(in) :: ssacld(:,:,:)     
    real(8), intent(in) :: asmcld(:,:,:)     
    real(8), intent(in) :: fsfcld(:,:,:)     
    real(8), intent(in) :: cicewp(:,:)       
    real(8), intent(in) :: cliqwp(:,:)       
    real(8), intent(in) :: reice(:,:)        
    real(8), intent(in) :: reliq(:,:)        
    real(8), intent(in) :: tauaer(:,:,:)     
    real(8), intent(in) :: ssaaer(:,:,:)     
    real(8), intent(in) :: asmaer(:,:,:)     
    real(8), intent(in) :: ecaer(:,:,:)      
    integer(4), intent(in) :: start_band     
    integer(4), intent(in) :: end_band       
    real(8), intent(inout) :: swuflx(:,:)      
    real(8), intent(inout) :: swdflx(:,:)      
    real(8), intent(inout) :: swhr(:,:)        
    real(8), intent(inout) :: swuflxc(:,:)     
    real(8), intent(inout) :: swdflxc(:,:)     
    real(8), intent(inout) :: swhrc(:,:)       
    real(8), intent(inout) :: wavenumber_range(:)

    integer icld_inout, iaer_inout
    
    icld_inout = icld
    iaer_inout = iaer
    call rrtmg_sw(                                      &
        ncol, nlay, icld_inout, iaer_inout,             &
        play, plev, tlay, tlev, tsfc,                   &
        h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr,   &
        asdir, asdif, aldir, aldif,                     &
        coszen, adjes, dyofyr, scon, isolvar,           &
        inflgsw, iceflgsw, liqflgsw, cldfr,             &
        taucld, ssacld, asmcld, fsfcld,                 &
        cicewp, cliqwp, reice, reliq,                   &
        tauaer, ssaaer, asmaer, ecaer,                  &
        start_band, end_band,                           &
        swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc,  &
        wavenumber_range,                               &
        bndsolvar, indsolvar, solcycfrac )

end subroutine rrtmg_sw_wrapper
