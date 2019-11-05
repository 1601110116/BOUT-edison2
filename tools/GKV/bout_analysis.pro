FUNCTION BOUT_test
gyro_data = xi_data()
p_n = gyro_data.p -> fft('z')
p_nsq = p_n -> abssq()
p_nsq -> get, axis='z', irange=n_range, gridrange=n_vrange
p_arr = objarr(n_range[1]+1)
p_arr2 = objarr(n_range[1]+1)
p_arr3 = objarr(n_range[1]+1)

dn = n_vrange[1]-n_vrange[0]
for i = 0,n_range[1] do begin
p_arr[i] = p_nsq -> slice(n=i*dn)
p_arr2[i] = p_arr[i] -> int('x')
p_arr3[i] = p_arr2[i] -> int('y')
endfor


RETURN, p_arr3
END
