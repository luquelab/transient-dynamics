def model_fun(t,y,beta_max,p,k,T,T_opt,gamma):
	w_daysy = y[0]
	b_daysy = y[1]

	w_daysy_dot = beta_max*p*w_daysy -beta_max*b_daysy*w_daysy -beta_max*w_daysy*w_daysy -p*k*(T-T_opt)**2*w_daysy -k*(T-T_opt)**2*b_daysy*w_daysy -k*(T-T_opt)**2*w_daysy*w_daysy -gamma*w_daysy 
	b_daysy_dot = beta_max*p*b_daysy -beta_max*w_daysy*b_daysy -beta_max*b_daysy*b_daysy -p*k*(T-T_opt)**2*b_daysy -k*(T-T_opt)**2*w_daysy*b_daysy -k*(T-T_opt)**2*b_daysy*b_daysy -gamma*b_daysy 
	
	y_dot = [w_daysy_dot ,b_daysy_dot ]
	return y_dot