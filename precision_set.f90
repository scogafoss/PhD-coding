module precision_set
	implicit none
	integer, parameter :: dp = selected_real_kind(8)
	integer, parameter :: re = selected_real_kind(4)
	integer, parameter :: qp = selected_real_kind(16)
	real(re),parameter :: pi_re = 4 * atan (1.0_re)
	real(dp),parameter :: pi_dp = 4 * atan (1.0_dp)
	real(qp),parameter :: pi_qp = 4 * atan (1.0_qp)
end module precision_set
