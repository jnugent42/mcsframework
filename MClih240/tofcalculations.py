//My formula TOF01 - good 
   dEdx = 3.00;
   double a = -1 + ((pow(c,2)*pow(t_initial,2))/pow(deltas,2));
   double b = -(((c*t_initial)*dEdx*(pow(5.287,2)-pow(12.92,2)+(2*12.92*16.95)-(2*5.287*16.95)))/pow(deltas,2));
   double cp = ((pow(pow(5.287,2)-pow(12.92,2)+(2*12.92*16.95)-(2*5.287*16.95),2)*pow(dEdx,2))/(4*pow(deltas,2)))-(pow(105.65,2));

   std::cout << "a " << a << std::endl;
   std::cout << "b " << b << std::endl;
   std::cout << "cp " << cp << std::endl;

   double holder = pow(b,2)-4*a*cp;
   std::cout << "holder " << holder << std::endl;
   double pcr = (-b+sqrt(pow(b,2)-(4*a*cp)))/(2*a);
   std::cout << "pcr " << pcr << std::endl;

//My formula TOF12 -  
   dEdx = -3.53;
   a = -1 + ((pow(c,2)*pow(t_12,2))/pow(21.13-12.92,2));
   b = -(((c*t_12)*dEdx*(pow(12.92,2)-pow(21.13,2)+(2*21.13*16.95)-(2*12.92*16.95)))/pow(21.13-12.92,2));
   cp = ((pow(pow(12.92,2)-pow(21.13,2)+(2*21.13*16.95)-(2*12.92*16.95),2)*pow(dEdx,2))/(4*pow(21.13-12.92,2)))-(pow(105.65,2));

   std::cout << "a " << a << std::endl;
   std::cout << "b " << b << std::endl;
   std::cout << "cp " << cp << std::endl;

   double pcr2 = (-b+sqrt(pow(b,2)-(4*a*cp)))/(2*a);
   std::cout << "pcr2 " << pcr2 << std::endl;a

// delta - from Bethe-Bloch (units are cm)
   double z = 3.25;
   double rho = 0.694;
   dEdx = BetheBloch(pz, 36.6e-6);
   std::cout << "dEdx " << dEdx << std::endl;
   double Eloss = dEdx * z * rho;
   double E = sqrt(pow(pz,2)+pow(105.65,2));
   E -= Eloss;
   double ploss = sqrt(pow(E,2) - pow(105.65,2));
   E += Eloss;
   double delta = pz - ploss;
   /*
   std::cout << "delta " << delta << std::endl;
   std::cout << "E " << E << std::endl;
   */

   // time0
   double t1 = s1*sqrt(pow(pz+delta,2)+pow(105.65,2))/(pz+delta);
   double t2 = s2*sqrt(pow(pz-delta,2)+pow(105.65,2))/(pz-delta);
   double t0 = t1 + t2;
   /*
   std::cout << "t2 " << t2 << std::endl;
   std::cout << "t1 " << t1 << std::endl;
   std::cout << "t0 " << t0 << std::endl;
   */

   // delta_t
   double term3i = 0.5* (s1 + s2)*pow(105.65,2)*pow(delta,2)/(pow(pz,2)*pow(E,2));
   double term3ii = pz/E + 2*E/pz;
   double term2 = -(s1 - s2) * delta * pow(105.65,2)/ (pow(pz,2)*E);
   double dtime = term2 + term3i*term3ii;
   /*
   std::cout << "term3i " << term3i << std::endl;
   std::cout << "term3ii " << term3ii << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "dtime " << dtime << std::endl;
   */

// Corrected P
   double time = t0 + dtime;
   p_corrected = pz - ((1/(s1-s2))*E*pow(pz,2)*dtime)/pow(105.65,2);
   std::cout << "pz " << pz << std::endl;
   pz = p_corrected;
   t_path = t0+dtime;
   std::cout << "p_corrected " << p_corrected << std::endl;
   /*
   std::cout << "p_corrected " << p_corrected << std::endl;
   std::cout << "counter " << counter << std::endl;
   std::cout << "iteration " << i << std::endl;
   std::cout << "t_path " << t_path << std::endl;
   */

// Energy loss number 
   double z = 3.25;
   double rho = 0.694;
   dEdx = BetheBloch(pz12,36.5e-6);
   double Eloss = dEdx * z * rho;
   std::cout << "Eloss " << Eloss << std::endl;
   z = 113;
   rho = 1.785e-4;
   dEdx = BetheBloch(pz12,41.8e-6);
   double ElossHe = dEdx * z * rho;
   std::cout << "ElossHe " << ElossHe << std::endl;
   z = 0.016;
   rho = 2.7;
   dEdx = BetheBloch(pz12,166.0e-6);
   double ElossAl = dEdx * z * rho;
   std::cout << "ElossAl " << ElossAl << std::endl;
   z = 0.74;
   rho = 1.04;
   dEdx = BetheBloch(pz12,68.7e-6);
   double ElossC8H8 = dEdx * z * rho;
   std::cout << "ElossC8H8 " << ElossC8H8 << std::endl;
   z = 15.24;
   rho = 1.032;
   dEdx = BetheBloch(pz12,68.7e-6);
   double ElossTOF = dEdx * z * rho;
   std::cout << "ElossTOF " << ElossTOF << std::endl;
   z = 4.6;
   rho = 0.2975;
   dEdx = BetheBloch(pz12,139e-6);
   double Elossckov = dEdx * z * rho;
   std::cout << "Elossckov " << Elossckov << std::endl;
   z = 764.2;
   rho = 1.23e-3;
   dEdx = BetheBloch(pz12,85e-6);
   double Elossair = dEdx * z * rho;
   std::cout << "Elossair " << Elossair << std::endl;
   std::cout << "E total " << ElossHe/2 +ElossAl/2 + ElossC8H8/2 + Eloss/2 + 1.5 + ElossTOF/2 + Elossckov<< std::endl;
   E -= ElossHe/2 +ElossAl/2 + ElossC8H8/2 + Eloss/2 + 1.5 + ElossTOF/2 + Elossckov;

std::cout << "t_12 " << t_12 << std::endl;
   // Paul Note formula 2nd order - OK
   dEdx = 3.53;
   double s_1 = 12.929-16.952;
   double s_2 = 21.139-16.952;
   double fterm = sqrt(1+((pow(dEdx,2)*pow(s_1+s_2,2))/(4*pow(105.65,2)*(pow(c*t_12/(8.21),2)-1))));
   double sterm = (dEdx*c*t_12*(s_1+s_2))/(2*8.21*(pow(c*(t_12)/(8.21),2)-1));
   double p12Paul2ndfor = pz12*(fterm)-sterm;

   std::cout << "fterm " << fterm << std::endl;
   std::cout << "sterm " << sterm << std::endl;

// Paul TOF12 long form 6th order formula - 
   dEdx = 3.53;
   term1 = pow(dEdx,4)*pow((pow(12.92,2)-pow(21.13,2)+(2*21.13*16.95)-(2*12.92*16.95)),2)/(4*pow(105,4));
   term2 = c*t_initial*pow(dEdx,3)*(pow(12.92,2)-pow(21.13,2)+(2*21.13*16.95)-(2*12.92*16.95))/(pow(105,3));
   term3 = -pow(dEdx,2)*pow(21.13-12.92,2)/(pow(105,2));
   term4 = pow(dEdx,2)*(pow(c,2)*pow(t_initial,2)-pow(21.13-12.92,2))/(pow(105,2));

   std::cout << "term1 " << term1 << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "term3 " << term3 << std::endl;
   std::cout << "term4 " << term4 << std::endl;

   x = 105/pz12;
   std::cout << "x " << x << std::endl;
   MyFunction1D myf3;
   myf3.a = term1;
   myf3.b = term2;
   myf3.cp = term3;
   myf3.d = term4;
   ROOT::Math::GradFunctor1D  f3(myf3);
   ROOT::Math::RootFinder rfn3(ROOT::Math::RootFinder::kGSL_NEWTON);
   rfn3.SetFunction(f3, x);
   rfn3.Solve();
   double p12longPaul6thfor = 105.5/rfn3.Root();
   std::cout << "p12longPaul6thfor " << p12longPaul6thfor << std::endl;

cor_mom->Fill(pcr);
double res12p2 = p12Paul2ndfor - MCTruth_pz_mid;
TOF12longPaul6thforvsMCTruth->Fill(p12longPaul6thfor,MCTruth_pz_mid);
   TOF01fordownvsMCTruth->Fill(pcr,MCTruth_pz_mid);
      TOF12fordownvsMCTruth->Fill(pcr2,MCTruth_pz_mid);
      double res01j = pcr - MCTruth_pz_mid;
      TOF12Paul2ndforvsMCTruth->Fill(p12Paul2ndfor,MCTruth_pz_mid);
      residual12p2->Fill(res12p2);
double res01p6 = p01Paul6thfor - MCTruth_pz_mid;
double rescobb = p_corrected-MCTruth_pz_mid;

// Paul TOF01 6th order formula 
   double dEdx = 3;
   double term1 = pow(dEdx,4)*pow((pow(5.287,2)-pow(12.92,2)+(2*12.92*16.95)-(2*5.287*16.95)),2)/(4*pow(105,4));
   double term2 = c*t_initial*pow(dEdx,3)*(pow(5.287,2)-pow(12.92,2)+(2*12.92*16.95)-(2*5.287*16.95))/(pow(105,3));
   double term3 = -pow(dEdx,2)*pow(12.92-5.287,2)/(pow(105,2));
   double term4 = pow(dEdx,2)*(pow(c,2)*pow(t_initial,2)-pow(12.92-5.287,2))/(pow(105,2));

   std::cout << "term1 " << term1 << std::endl;
   std::cout << "term2 " << term2 << std::endl;
   std::cout << "term3 " << term3 << std::endl;
   std::cout << "term4 " << term4 << std::endl;

   double x = 105/pz01;

   MyFunction1D myf1;
   myf1.a = term1;
   myf1.b = term2;
   myf1.cp = term3;
   myf1.d = term4;
   std::cout << "myf1.operator()(pz) " << myf1.operator()(x) << std::endl;
   std::cout << "myf1.Derivative(pz) " << myf1.Derivative(x) << std::endl;
   ROOT::Math::GradFunctor1D  f1(myf1);
   ROOT::Math::RootFinder rfn(ROOT::Math::RootFinder::kGSL_NEWTON);
   rfn.SetFunction(f1, x);
   if (isfinite(myf1.operator()(pz)) && isfinite(myf1.Derivative(x))) {
   rfn.Solve();
   }
   double p01Paul6thfor = 105.5/rfn.Root();
   std::cout << "p01Paul6thfor " << p01Paul6thfor << std::endl;

