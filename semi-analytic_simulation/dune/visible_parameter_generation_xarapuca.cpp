//
// root macro to create semi-analytic model parameters for visible photons
// DUNE, X-Arapucas

struct acc{
	// ax,ay,az = centre of rectangle; w = y dimension; h = z dimension
	double ax, ay, az, w, h; 
};

// functions
int VisHits(const int Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint);
Double_t GaisserHillas(double x,double *par);
double omega(const double &a, const double &b, const double &d);
double solid(const acc& out, const TVector3 &v);
double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate );

double GHslope(double theta, double p0, double p1) {
  double y = p1*theta + p0;
  return y;
}

// required constants:
// optical detector information
const double y_dimension_detector = 9.3;	// cm
const double z_dimension_detector = 46.8;	// cm

// cathode foils information - DUNE
// tpc drift length
const double plane_depth = 363.38405;	// cm 
// size of cathode covered by foils
const double y_dimension_foils = 1204.7255 + 5.466;	// cm		// 2 panels y height 602.36275 with 5.466cm gaps between
const double z_dimension_foils = 1359.144 + 35.196;	// cm		// 6 panels of z width 226.524cm with 5.866cm gaps between them
// centre coordinates of foils
const double x_foils = 363.38405; const double y_foils = 0; const double z_foils = 696.294;	// cm

// LAr properties
const double L_abs = 2000;
const double scint_light_yield = 24000; // photons/MeV

const double pi = 3.1416;

// VUV Gaisser-Hillas parameters

// DUNE flat surface (for cathode)
// DUNE - RS99cm
const double fGHVUVPars[4][9] = { {1.23286, 1.20491, 1.1597, 1.08148, 0.986607, 0.868075, 0.725637, 0.633551, 0.469717},
  								{150.325, 150.109, 156.294, 169.271, 179.642, 212.9, 232.173, 226.557, 215.166},
  								{19.0558, 18.8478, 18.8757, 18.7149, 50.8662, 49.2247, 57.6717, 157.92, 172.494},
  								{-3000, -3000, -3000, -3000, -1000, -1000, -1000, -100, -100} };
vector<double> angulo = {5, 15, 25, 35, 45, 55, 65, 75, 85};
vector<double> slopes1 = {-9.63855e-05, -6.82604e-05, -9.63478e-05, -0.000121181, -0.000126611, -0.000115481, -8.61492e-05, -0.000112594, -7.80935e-05};
vector<double> slopes2 = {-0.0662469, -0.0504497, -0.0596321, -0.0418021, -0.0342462, -0.0531668, -0.0522639, -0.0578887, -0.0591081};
vector<double> slopes3 = {-0.00593207, -0.00672713, -0.0020843, -0.00216374, 0.00901291, 0.00385402, 0.0066081, 0.0341547, 0.0446519};


// fit settings
// angular bin size, deg								  
const int number_angle_bins = 9;
const double delta_angulo = 90 / number_angle_bins;

// range and step of profiling
const double d_min = 0;
const double d_max = 365.;
const double step_d = 1;

// optical detector list
const string pmtlistname = "../../optical_detectors/opch_info_dune_1x2x6.txt"; // DUNE geometry

// input filename
// most recent file located on /dune/data/users/pgreen/optical_parameterisations/
string inputfilename = "../../data/dune_1x2x6_withWires_withFoils.root";

// chosen subset in r (radial distance from centre) -- select range around central value to capture the line of points to get full distribution
// check number of files selected makes sense / there are not holes in distribution
const double r_range[2] = {0, 75};
//const double r_range[2] = {140, 160};
//const double r_range[2] = {240, 275};
//const double r_range[2] = {340, 390};
//const double r_range[2] = {485, 505};
//const double r_range[2] = {575, 625};
//const double r_range[2] = {670, 720};
//const double r_range[2] = {780, 805};
//const double r_range[2] = {825, 1000};

const double min_number_entries = 250;	 // matching VUV

// function to calculate corrections for borders
void visible_parameter_generation_xarapuca(){

	gRandom->SetSeed(0);

	// load optical detector positions from text file
	vector<vector<double>> optical_detector_positions;
	cout << "Loading Photon Detector positions...\n";
    ifstream myfile;
    myfile.open(pmtlistname); 
    if(myfile.is_open()) {
    	cout << "File opened successfully" << endl; 
	    while(!myfile.eof()) {
	        double num_pmt, x_pmt, y_pmt, z_pmt;
	        if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt) {     
	            	vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
	            	optical_detector_positions.push_back(line_data);
	        }
	        else break;
	    }
	    myfile.close();
	    cout << "Photon detector positions loaded.\n\n";
	}
	else cout << "Optical detector list not found." << endl;

	int numberPMTs = optical_detector_positions.size();
	
	cout<<"----> number optical detectors: " << numberPMTs << endl;
  	
    // vectors to store calulated values
    vector<double> v_distance; v_distance.reserve(1e6);
    vector<double> v_hits_sim; v_hits_sim.reserve(1e6);
    vector<double> v_hits_geo; v_hits_geo.reserve(1e6);
    vector<double> v_offset_angle; v_offset_angle.reserve(1e6);
    vector<double> v_r; v_r.reserve(1e6);
   
    // open file
    TFile* f = new TFile(inputfilename.c_str());
	TTree *tree = (TTree *)f->Get("myTree");
	const Int_t kMaxDevices = 600;
	int VUV_hits[kMaxDevices];
	int Vis_hits[kMaxDevices];
	double X, Y, Z;
	int genPhotons, numberDevices;
	tree->SetBranchAddress("numberDevices",&numberDevices);
	tree->SetBranchAddress("X", &X);
	tree->SetBranchAddress("Y", &Y);
	tree->SetBranchAddress("Z", &Z);
	tree->SetBranchAddress("VUV_hits",VUV_hits);
	tree->SetBranchAddress("Vis_hits",Vis_hits);
	tree->SetBranchAddress("genPhotons", &genPhotons);
	    
	// loop through TTree
	cout << tree->GetEntries() << endl;
	for(int n=0; n < tree->GetEntries(); n++) {
		tree->GetEntry(n);
    	double posSource[3] = {X, Y, Z};

    	// check if within range of chosen parameter set
	    double R = sqrt( pow(posSource[1],2) + pow(posSource[2] - z_foils, 2) );
	    if (R < r_range[0] || R > r_range[1]) continue;

		cout << "Scintillation point --> (x, y, z) = " << posSource[0] << "  " << posSource[1] << "  " << posSource[2] << ", R = " << R << endl;

		// ****
		// Calculate amount of light: 
		// ****
		
		// loop through optical channels
		for (int nPMT = 0; nPMT < numberPMTs; nPMT++) {
			
			// get Scintpoint
			TVector3 ScintPoint(posSource[0],posSource[1],posSource[2]);
			// get OpDetPoint
			int pmt_index = optical_detector_positions.at(nPMT).at(0);
			TVector3 OpDetPoint(optical_detector_positions.at(nPMT).at(1),optical_detector_positions.at(nPMT).at(2),optical_detector_positions.at(nPMT).at(3));
			
			// get hotspot
			TVector3 hotspot(plane_depth, ScintPoint[1], ScintPoint[2]);

			// calculate vis distance and angle [between hotspot and opdetpoint]
			double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
			double cosine = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
	  		double theta = acos(cosine)*180./pi;

	  		// distance to cathode
	  		double x_distance = posSource[0];

			// solid angle prediction
			double nPhotons_solid = VisHits(genPhotons,ScintPoint,OpDetPoint);

			if (Vis_hits[pmt_index] < min_number_entries) {
				cout << "Point Removed: Low Stats, N = " << Vis_hits[pmt_index] << endl;
				continue;
			}

			// save values
			v_distance.push_back(plane_depth - x_distance);			// distance, x position
			v_offset_angle.push_back(theta);			// offset angle
			v_hits_sim.push_back(Vis_hits[pmt_index]);	// full simulation visible hits
			v_hits_geo.push_back(nPhotons_solid);		// solid angle visible hits
			v_r.push_back(R);							// radial distance from centre
			
		} // end loop through optical channels
		
	} // end of loop over points
	
	// create profiles, one for each angular bin	
	TProfile* pdiff_d[number_angle_bins];
	for(int loop = 0; loop < number_angle_bins; loop++) {
		double theta_min = loop*delta_angulo;
		double theta_max = theta_min + delta_angulo;
		pdiff_d[loop]  = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
	}

	// populate profiles
	vector<int> nEntries(9,0);
	for(int i=0; i < v_distance.size(); i++) {
	    double costheta = cos(pi*v_offset_angle.at(i)/180.);
	    // which angular bin
	    int j = int(v_offset_angle.at(i)/delta_angulo);
	    pdiff_d[j]->Fill(v_distance.at(i), v_hits_sim.at(i)/v_hits_geo.at(i)*costheta);
	    nEntries[j]++;
 	}

 	// find mean value of r
	double sum = std::accumulate(std::begin(v_r), std::end(v_r), 0.0);
	double r =  sum / v_r.size();
	
 	// create plot
 	TCanvas *c1 = new TCanvas("c1","",200,200,1000,1000);
	c1->SetBottomMargin(0.105);
	c1->SetLeftMargin(0.105);

    auto leg = new TLegend(0.55, 0.4 + 0.0544*1 + 0.05, 0.89, 0.89 + 0.05, NULL, "brNDC");
    leg->SetHeader("");
  	leg->SetBorderSize(0);
  	leg->SetFillStyle(0);

  	char label[number_angle_bins][20];

    // plot to format axis
    double x_0[2] = {0, 365};
  	double y_0[2] = {0, 4.6};
  	TGraph* gg = new TGraph(2,x_0,y_0);
  	gg->GetXaxis()->SetTitle("d_{c} [cm]");
	gg->GetYaxis()->SetTitle("N_{Geant4} / N_{#Omega} / cos(#theta_{c})");
	//gg->GetYaxis()->SetLabelSize(0.05);
	gg->GetYaxis()->SetTitleSize(0.05);
	gg->GetYaxis()->SetTitleOffset(1.);
	gg->GetYaxis()->SetRangeUser(0,4.6);
	
	//gg->GetXaxis()->SetLabelSize(0.05);
	gg->GetXaxis()->SetTitleSize(0.05);
	gg->GetXaxis()->SetTitleOffset(1.);
	gg->GetXaxis()->SetRangeUser(0,365);

	char title[100];
	sprintf(title,"DUNE XArapucas");
	gg->SetTitle(title);
	gg->Draw("ap");

 	// loop through profiles creating plots
 	for(int loop = number_angle_bins - 1; loop >= 0; loop--) {		
 		if (nEntries[loop] == 0) continue; 		

 		// profile
 		if (loop == 4) {
 			pdiff_d[loop]->SetLineColor(kOrange+7);
	    	pdiff_d[loop]->SetMarkerColor(kOrange+7);
 		}
 		else {
 			pdiff_d[loop]->SetLineColor(1 + loop);
	    	pdiff_d[loop]->SetMarkerColor(1 + loop);
	    }

	    pdiff_d[loop]->SetMarkerStyle(20 + loop);
	    pdiff_d[loop]->SetMarkerSize(0.9);
	    pdiff_d[loop]->Draw("E1 p same");	    
 	}

 	
 	for(int loop = 0; loop < number_angle_bins; loop++) {		
 		if (nEntries[loop] == 0) continue; 

 		// legend
	    int a_min = loop*delta_angulo;
    	int a_max = (loop+1)*delta_angulo;
	    sprintf(label[loop],"#theta_{c} #in [%i, %i] deg",a_min, a_max);
	    leg->AddEntry(pdiff_d[loop],label[loop],"p");	    

	}

 	leg->Draw("same");


 	// extract parameters from profiles and output
 	vector<vector<double>> pars;
 	vector<vector<double>> pars_errs;
 	vector<double> x_positions;
 	for(int loop = 0; loop < number_angle_bins; loop++) {	
 		if (nEntries[loop] == 0){
 			// skip empty angle bins, not all large angles exist due to geometry
 		 	// parameters set to parameter for previous bin
 		 	pars.push_back(pars[loop-1]);
 		 	continue;
 		}

 		vector<double> pars_loop;
 		vector<double> pars_errs_loop;
 		int n_bins = (d_max-d_min)/step_d;
 		for (int i = 0; i < n_bins; i++){
 			// skip empty bins
    		if(!(pdiff_d[loop]->GetBinEntries(i) > 0)) continue;
    		// save parameters
    		pars_loop.push_back(pdiff_d[loop]->GetBinContent(i));
    		pars_errs_loop.push_back(pdiff_d[loop]->GetBinError(i));
    		// save x positions on first loop
    		if (loop == 0) {
    			x_positions.push_back((i-1) * step_d);
    		}    	 	
    	}
    	pars.push_back(pars_loop);
    	pars_errs.push_back(pars_errs_loop);
 		
 	}  

 	cout << "R = " << r << endl;
 	cout << "X values = ";
 	for (int i = 0; i < x_positions.size(); i++) {
 		if (i == x_positions.size() - 1){
 			cout << x_positions[i] << endl;
 		}
 		else cout << x_positions[i] << ", ";
 	}


 	cout << "Correction Parameters: " << endl;
 	for (int angle_bin = 0; angle_bin < pars.size(); angle_bin++){
 		for (int x_bin = 0; x_bin < pars[0].size(); x_bin++){
 			if (x_bin == pars[0].size() - 1){
 				cout << pars[angle_bin][x_bin] << endl;
 			}
 			else cout << pars[angle_bin][x_bin] << ", ";
 		}
 	}

 	cout << "Correction Parameters Uncertainties: " << endl;
 	for (int angle_bin = 0; angle_bin < pars_errs.size(); angle_bin++){
 		for (int x_bin = 0; x_bin < pars_errs[0].size(); x_bin++){
 			if (x_bin == pars_errs[0].size() - 1){
 				cout << pars_errs[angle_bin][x_bin] << endl;
 			}
 			else cout << pars_errs[angle_bin][x_bin] << ", ";
 		}
 	}    
}


// visible number of hits calculation function
// used for parameter generation - does not apply any correction
// Visible hits calculation
int VisHits(const int Nphotons_created, const TVector3 &ScintPoint, const TVector3 &OpDetPoint) {
 
	// 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:

	// set cathode plane struct for solid angle function
	acc cathode_plane; 
	cathode_plane.ax = x_foils; cathode_plane.ay = y_foils; cathode_plane.az = z_foils;   		// centre coordinates of cathode plane
	cathode_plane.w = y_dimension_foils; cathode_plane.h = z_dimension_foils;                     // width and height in cm

	// get scintpoint coords relative to centre of cathode plane
	TVector3 cathodeCentrePoint(x_foils,y_foils,z_foils);
	TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint; 

	// calculate solid angle of cathode from the scintillation point
	double solid_angle_cathode = solid(cathode_plane, ScintPoint_relative);

	// calculate distance and angle between ScintPoint and hotspot
	// vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
	double distance_cathode = std::abs(x_foils - ScintPoint[0]); 

	double cosine_cathode = 1;
	double theta_cathode = 0.0;

	// calculate hits on cathode plane via geometric acceptance
	double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4.*pi)) * Nphotons_created;

	// apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
	// offset angle bin
	int j = (theta_cathode/delta_angulo);  
	// correction
	double pars_ini[4] = {fGHVUVPars[0][j], fGHVUVPars[1][j], fGHVUVPars[2][j], fGHVUVPars[3][j]};

	// gh border
	double r_distance = sqrt( pow(ScintPoint[1] - y_foils, 2) + pow(ScintPoint[2] - z_foils, 2)); 

	double s1 = interpolate( angulo, slopes1, theta_cathode, true);
	double s2 = interpolate( angulo, slopes2, theta_cathode, true);
	double s3 = interpolate( angulo, slopes3, theta_cathode, true);

	pars_ini[0] = pars_ini[0] + s1 * r_distance;
	pars_ini[1] = pars_ini[1] + s2 * r_distance;
	pars_ini[2] = pars_ini[2] + s3 * r_distance;
	pars_ini[3] = pars_ini[3];

	double GH_correction = GaisserHillas(distance_cathode, pars_ini);

	double cathode_hits_rec = GH_correction*cathode_hits_geo/cosine_cathode;


	// 2). calculate number of these hits which reach the optical channel from the hotspot via solid angle 

	// calculate hotspot location  
	TVector3 v_to_wall(x_foils - ScintPoint[0],0,0);        
	TVector3 hotspot = ScintPoint + v_to_wall;

	// distance to hotspot
	double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
	// distance from hotspot to arapuca
	double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
	// angle between hotspot and arapuca
	double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
	double theta_vis = acos(cosine_vis)*180./pi;

	// solid angle :
	double solid_angle_detector = 0;
	// rectangular aperture  
	// set Arapuca geometry struct for solid angle function
	acc detPoint; 
	detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
	detPoint.w = y_dimension_detector; detPoint.h = z_dimension_detector; // width and height in cm of arapuca active window

	// get hotspot coordinates relative to detpoint
	TVector3 emission_relative = hotspot - OpDetPoint;

	// calculate solid angle of optical channel
	solid_angle_detector = solid(detPoint, emission_relative);

	// calculate number of hits via geometeric acceptance  
    double hits_geo = (solid_angle_detector / (2*pi)) * cathode_hits_rec;

	// round final result
	int hits_vis = std::round(hits_geo);

	return hits_vis;
}

// semi-analytic model utility functions
// gaisser-hillas function definition
Double_t GaisserHillas(double x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-x)/par[2]);
  
  return ( Normalization*Term*Exponential);
}

// solid angle of rectanglular aperture calculation functions
double omega(const double &a, const double &b, const double &d){

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double solid(const acc& out, const TVector3 &v){

  //v is the position of the track segment with respect to 
  //the center position of the arapuca window 
 
  // arapuca plane fixed in x direction	

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

double interpolate( const std::vector<double> &xData, const std::vector<double> &yData, double x, bool extrapolate ) {
  int size = xData.size();
  int i = 0;                                          // find left end of interval for interpolation
  if ( x >= xData[size - 2] )                         // special case: beyond right end
    {
      i = size - 2;
    }
  else
    {
      while ( x > xData[i+1] ) i++;
    }
  double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1]; // points on either side (unless beyond ends)
  if ( !extrapolate )                                                    // if beyond ends of array and not extrapolating
    {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
    }
  double dydx = ( yR - yL ) / ( xR - xL );            // gradient
  return yL + dydx * ( x - xL );                      // linear interpolation
}