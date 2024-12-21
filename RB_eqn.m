%% 6DOF rigid body dynamics equation in DQ form 
function[result] = RB_eqn (w_dq_b_i_b,J_dq,F)
A = cross((w_dq_b_i_b),J_dq*(w_dq_b_i_b));
result = inv(J_dq)*(F - A);
end
