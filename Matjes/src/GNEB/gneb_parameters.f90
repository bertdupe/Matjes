module m_gneb_parameters

character(len=35) :: momfile_i                             !< File name for initial magn. configuration
character(len=35) :: momfile_f                             !< File name for final magn. configuration
character(len=35) :: restartfile_if                        !< Name of the file with initial and final states needed to restart energy minimization calculations
character(len=35) :: restartfile_path                      !< Name of the file with the path needed to restart GNEB calculations
real(kind=8) :: ind_tol
real(kind=8) :: amp_rnd                                   !< Amplitude of random perturbation of the components of magnetic moments
real(kind=8) :: amp_rnd_path
integer :: minalgo                                           !< Minimization algorithm (1=VPO,...)
integer :: minitrmax                                         !< Maximum number of iterations in energy minimization procedure
real(kind=8) :: minftol                                     !< Convergence criterion in mRy
integer :: mintraj_step                                      !< Save configuration every 'mintraj_step' step during energy minimization
real(kind=8) :: vpodt                                       !< VPO timestep
real(kind=8) :: vpomass                                     !< VPO mass
integer :: initpath                                          !< Path initialization method (1=geodesic, 2=read from file,...)
real(kind=8) :: spring                                      !< Magnitude of the spring constant in GNEB calculations
integer :: mepitrmax                                         !< Maximum number of iterations in MEP finding procedure
real(kind=8) :: mepftol                                     !< Convergence criterion in the GNEB method, in mRy
real(kind=8) :: mepftol_ci                                  !< Convergence criterion in the CI-GNEB method, in mRy
integer :: meptraj_step                                      !< Save configuration every 'meptraj_step' step during MEP finding procedure
character(len=1) :: do_gneb                                  !< Do GNEB calculations (Y/N)
character(len=1) :: do_gneb_ci                               !< Do CI-GNEB calculations (Y/N)
character(len=1) :: do_norm_rx                               !< Normalize reaction coordinate (Y/N)
character(len=1) :: en_zero                                  !< Level of zero energy. 'I' - initial state; 'F' - final state; 'N' - 0.0
integer :: sample_num                                         !< Number of samples in the interpolated curve
integer :: nim

public

contains

subroutine read_gneb_parameters()
use m_io_files_utils
use m_io_utils
implicit none
integer :: io
character(len=30) :: fname

fname='input'
!    character(len=1) :: OPT_flag_str, ip_adapt_flag_str, OPT_printcores_flag_str

io=open_file_read(fname)

call get_parameter(io,fname,'momfile_i',momfile_i)
call get_parameter(io,fname,'momfile_f',momfile_f)
call get_parameter(io,fname,'restartfile_if',restartfile_if)
call get_parameter(io,fname,'restartfile_path',restartfile_path)

call get_parameter(io,fname,'spring',spring)
call get_parameter(io,fname,'initpath',initpath)
call get_parameter(io,fname,'amp_rnd',amp_rnd)
call get_parameter(io,fname,'amp_rnd_path',amp_rnd_path)
call get_parameter(io,fname,'min_itrmax',minitrmax)
call get_parameter(io,fname,'mintraj_step',mintraj_step)
call get_parameter(io,fname,'min_ftol',minftol)
call get_parameter(io,fname,'mep_itrmax',mepitrmax)
call get_parameter(io,fname,'meptraj_step',meptraj_step)
call get_parameter(io,fname,'mep_ftol',mepftol)
call get_parameter(io,fname,'mep_ftol_ci',mepftol_ci)
call get_parameter(io,fname,'io_do_gneb',do_gneb)
call get_parameter(io,fname,'do_gneb_ci',do_gneb_ci)
call get_parameter(io,fname,'do_norm_rx',do_norm_rx)
call get_parameter(io,fname,'en_zero',en_zero)
call get_parameter(io,fname,'vpo_dt',vpodt)
call get_parameter(io,fname,'vpo_mass',vpomass)
call get_parameter(io,fname,'sample_num',sample_num)
call get_parameter(io,fname,'nim',nim)

call close_file(fname,io)

end subroutine read_gneb_parameters

subroutine set_gneb_defaults()

implicit none
momfile_i     = 'momfile_i'
momfile_f     = 'momfile_f'
restartfile_if = 'nofile'
restartfile_path = 'nofile'
amp_rnd = 0.0
amp_rnd_path = 0.0

minftol = 0.000000001d0
mintraj_step = 100
vpodt = 0.01d0
vpomass = 1d0
minitrmax = 10000000

      !Parameters for GNEB calculations
initpath = 1
spring = 0.5d0
mepftol = 0.001d0
mepftol_ci = 0.00001d0
mepitrmax = 10000000
meptraj_step = 10
do_gneb = 'Y'
do_gneb_ci = 'N'
do_norm_rx = 'N'
en_zero = 'N'
nim = 10

    !Parameters for energy interpolation along the MEP
sample_num = 500
end subroutine set_gneb_defaults

end module m_gneb_parameters
