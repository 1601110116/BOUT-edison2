	if (MYPE==0) {
		output.write("output all the polygons: domain,solcore,pf1,pf2,core.\n");
		output.write("%5i\n",ndomain);
		for (int i=0; i<ndomain;i++) {
			output.write("%20.10e%20.10e\n",r_domain(i),z_domain(i));
		}
		output.write("%5i\n",nsolcore);
		for (int i=0; i<nsolcore;i++) {
			output.write("%20.10e%20.10e\n",r_solcore(i),z_solcore(i));
		}
		output.write("%5i\n",npf1);
		for (int i=0; i<npf1;i++) {
			output.write("%20.10e%20.10e\n",r_pf1(i),z_pf1(i));
		}
		output.write("%5i\n",npf2);
		for (int i=0; i<npf2;i++) {
			output.write("%20.10e%20.10e\n",r_pf2(i),z_pf2(i));
		}
		output.write("%5i\n",ncore);
		for (int i=0; i<ncore;i++) {
			output.write("%20.10e%20.10e\n",r_core(i),z_core(i));
		}
	}
