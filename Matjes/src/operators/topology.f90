module m_topology_commons
use m_derived_types, only : operator_real
! operator that returns the topological charge
type(operator_real),save :: charge
! operator that returns the vorticity
type(operator_real),save :: vorticity

end module m_topology_commons
