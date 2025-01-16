#define WATER             //!< new in v.3.0: If defined, an explicit water cycle is added, with an explicit belowground space. The horizontal resolution of the soil field is currently set by DCELL


int nblayers_soil; //!< Global variable: number of soil layers (for water module)
float *layer_depth(0); //!< Global vector: depth of each layer (m) !!!UPDATE

char buffer[256], inputfile[256], inputfile_daytimevar[256], inputfile_climate[256], inputfile_soil[256], outputinfo[256],inputfile_inventory[256], inputfile_pointcloud[256], *bufi(0), *bufi_daytimevar(0), *bufi_climate(0), *bufi_soil(0), *buf(0), *bufi_data(0), *bufi_pointcloud(0); //!< Global variable: character strings used to read file names, or other features
#ifdef WATER
char inputfile_SWC[256], *bufi_dataSWC(0);

int _SOIL_LAYER_WEIGHT; // !< User control: three ways of computing the fraction of transpiration supplied by each soil layer in individual tree total transpiration, and to weight the tree average water potential in the root zone. (0: root biomass only, 1: relative root-to soil conductance, 2: relative estimated maximal transpiration as in Duursma & Medlyn 2012).
int _WATER_RETENTION_CURVE; // !< User control: different water retention curves can be used. So far two are have been implemented: either brooks & Corey (0), either van Genuchten-Mualem (1). To each wtare retention cirve option is associated a different set of pedo-transfer functions, see Table 2 (texture-based Tomasella & Hodnett 1998


// soil parameters (Sat_SWC, Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this should be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils, and comments in the code). Alternatively, if no local field soil data is available, these soil parameters (Sat_SWC, Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (Sat_SWC, Res_SWC) should probably be provided in input, and the computation of those properties from the available local data made using a new function of RconTROLL, if unearthed.
// since soil layers (silt, clay, sand) are only needed locally, they are now coded as vectors
//float *proportion_Silt(0);             //!soil layer silt fraction
//float *proportion_Clay(0);             //!soil layer clay fraction
//float *proportion_Sand(0);             //!soil layer sand fraction
float *Sat_SWC(0);          //!< Global vector: soil layer saturated water content, in m3/m3 -- this is often assumed similar to porosity, even though it is usually 5-10% lower than total porosity due to entrapped or dissolved air -- see comment Table 1 in Marthews et al. 2014
float *Max_SWC(0);          //!< Global vector: soil layer maximum absolute water content (m^3)
float *FC_SWC(0);          //!< Global vector: soil layer maximum absolute water content (m^3)
float *Res_SWC(0);          //!< Global vector: soil layer residual water content (m^3/m^3)
float *Min_SWC(0);          //!< Global vector: soil layer minimum absolute water content (m^3)
float *Ksat(0);             //!< Global vector: soil layer saturated conductivity (m^3)
float *a_vgm(0);            //!< Global vector: parameter for the van Genuchten_mualem soil water retention curves
float *b_vgm(0);            //!< Global vector: parameter for the van Genuchten_mualem soil water retention curves
float *c_vgm(0);            //!< Global vector: parameter for the van Genuchten_mualem soil water retention curves
float *m_vgm(0);            //!< Global vector: parameter for the van Genuchten_mualem soil water retention curves
float *phi_e(0);            //!< Global vector: parameter for the Campbell-Mualem soil water retention curves (possible update: replace with a Genuchten parameter)
float *b(0);                //!< Global vector: parameter for the Campbell-Mualem soil water retention curves (possible update: replace with a Genuchten parameter)
float **SWC3D(0);           //!< Global 3D field: soil water content in each soil voxel (layer * DCELL)
float **soil_phi3D(0);      //!< Global 3D field: soil water potential (in MPa) in each soil voxel (layer * DCELL)
float **Ks(0);              //!< Global 3D field: soil hydraulic conductivity in each soil voxel (layer * DCELL)
float **KsPhi(0);           //!< Global vector: soil hydraulic conductivity * soil water potential for each soil voxel (layer * DCELL), useful to ease computation
float **LAI_DCELL(0);        //!< Global vector: total leaf area index (LAI), averaged per DCELL
float *LAI_young(0);        //!< Global vector: total young leaf area index (LAI), averaged across all sites
float *LAI_mature(0);        //!< Global vector: total mature leaf area index (LAI), averaged across all sites
float *LAI_old(0);        //!< Global vector: total old leaf area index (LAI), averaged across all sites
float *Canopy_height_DCELL(0); //!< Global vector: mean top canopy height, averaged per DCELL
int *HSum_DCELL(0);       //!< Global vector: number of sites covered by vegetation per DCELL, used to compute mean top canopy height per DCELL
float *TopWindSpeed_DCELL(0); //!< Global vector: wind speed at the top of the canopy, per DCELL
float *Interception(0);     //!< Global vector: water interception by the canopy, per DCELL !!!UNITS
float *Throughfall(0);      //!< Global vector: throughfall, per DCELL !!!UNITS
float *Runoff(0);           //!< Global vector: water run-off, per DCELL !!!UNITS
float *Leakage(0);          //!< Global vector: water leakage, per DCELL !!!UNITS
float *Evaporation(0);      //!< Global vector: water evaporation (physical process), per DCELL
float **Transpiration(0);   //!< Global vector: water uptake by trees, in each soil voxel (layer * DCELL)
float transpiration_1016; //variable used only to compare with eddy-flux tower data.

float abund_phi_root; //!< summary statistic output: abundance-weighted phi-root
float abund10_phi_root; //!< summary statistic output: abundance-weighted phi-root, including only tree with dbh>10cm
float agb_phi_root; //!< summary statistic output: agb-weighted phi-root
#endif

// !!!: suggestion, maybe redefine some quantities via a field (i.e. tree root biomass...) or locally where they are needed (especially if they need to be recomputed every timestep)
float t_root_depth;     //!< Tree rooting depth (m)
float t_phi_root;       //!< Soil water potential in the root zone (MPa)
vector<float> t_root_biomass;   //!< Tree root biomass in each soil layer, in g
vector<float> t_soil_layer_weight;  //!< Soil layer weight (to compute t_phi_root), if different from root biomass
float t_WSF;            //!< Tree water stress factor for stomatal conductance, unitless, between 0 and 1
float t_WSF_A;          //!< Tree water stress factor for photosynthetic capacities, unitless, between 0 and 1
float t_transpiration;  //!< Amount of water taken up from the soil and transpired at each timestep !!!UPDATE: which unit?
float t_g1_0;
float t_g1;

//! Function constructor Tree()
    Tree(){
        t_from_Data = 0;
        t_sp_lab = 0;
        t_age = 0;
        t_hurt = 0;
        t_NPP=t_GPP=t_Rday=t_Rnight=t_Rstem=0.0; // new v.2.2
        t_dbh = t_height = t_CR = t_CD= 0.0;
        t_CrownDisplacement = 0;
#ifdef MIP_Lichstein
        t_inInventory=0;
#endif
        
        if(_NDD){
            t_NDDfield.reserve(nbspp+1);
            for(int sp=0;sp<(nbspp+1);sp++) t_NDDfield.push_back(0.0);
        }
#if defined(WATER)
        t_transpiration = 0.0;
        t_root_biomass.reserve(nblayers_soil);
        for(int l=0;l<nblayers_soil;l++) {
            t_root_biomass.push_back(0.0);
        }
        t_soil_layer_weight.reserve(nblayers_soil);
        for(int l=0;l<nblayers_soil;l++){
            t_soil_layer_weight.push_back(0.0);
        }

#ifdef WATER
    void Water_availability();      //!< Computation of the tree water availability in the root zone
    //compute root depth, root biomass in each layer, soil water potential in the root zone, and water stress factor
    //see comments at Tree::Water_availability
    //void UpdateRootDistribution();  //compute root depth, root biomass in each layer, soil water potential in the root zone, and water stress factor
    void Water_uptake();            //!< Contribution of trees to the stand Transpiration field -- called by UpdateField
#endif


#ifdef WATER
        //*###########*/
        //*## water ##*/
        //*###########*/
        
        t_transpiration = 0.0;
        
        parameter_name = "transpiration";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_transpiration, 0.0f, 1.0f, 0.0f, quiet);

        parameter_name = "root_depth";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_root_depth, 0.0f, 100.0f, 2.0f, quiet);
        
        parameter_name = "phi_root";
        parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
        SetParameter(parameter_name, parameter_value, t_phi_root, -15.0f, 0.0f, 0.0f, quiet);
        
        for(int l=0;l<nblayers_soil;l++){
            char par_name[30];
            sprintf(par_name, "root_biomass%d", l);
            string par_name_s = par_name;
            parameter_name = par_name_s;
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_root_biomass[l], 0.0f, 100000.0f, 1.0f, quiet);
        }
        
        
        for(int l=0;l<nblayers_soil;l++){
            char par_name[30];
            sprintf(par_name, "soil_layer_weight%d", l); //
            string par_name_s = par_name;
            parameter_name = par_name_s;
            parameter_value = GetParameter(parameter_name, parameter_names, parameter_values);
            SetParameter(parameter_name, parameter_value, t_soil_layer_weight[l], 0.0f, 1.0f, 1.0f, quiet);

void Tree::Water_availability() {
    
    //if(t_LA > 0.0){       // TO BE CHECKED: what happens to a tree with leafarea=0 in the following timesteps ?
    //Tree root biomass
    float total_root_biomass=t_LA*t_LMA;
    //In the absence of a more explicit allocation scheme, we assume that the total root biomass=the total leaf biomass. This may be discussed and changed. Eg. in Schippers et al. 2015 Functional Plant Biology, they assume that root biomass is 1.4 times leaf biomass... More importantly, the tree-level carbon budget should be closed and balanced. Needs to be specified that this correponds to the total FINE root biomass. Note that similarly, in ED2 (Xu et al. 2016), a constant ratio of 1 is assumed between leaf biomass and fine root biomass (whereas a constant ratio of 0.25 is assumed between coarse root and aboveground stem biomass). But leaf and fine root phenology follow slightly different dynamics (see Xu et al. 2016 SI) so unclear how this is consistent...?
    
    //Tree root depth
    //t_root_depth=fminf(layer_depth[(nblayers_soil-1)],pow((0.001*total_root_biomass),0.8)*0.2173913);       //this is the tree root depth, computed according to Arora & Boer 2003 Earth interactions (see equ. 11 therein). Instead of assuming a fixed root exponential profile, that typically results in a fixed root depth from recruitment to death, hence overestimating root vertical growth and depth for small saplings, they derived a root distribution profile that depends on root biomass. The simulated root depth and distribution have to be checked. 0.2173913=3/b with b=13.8
    //t_root_depth=fminf(layer_depth[(nblayers_soil-1)],pow((0.001*total_root_biomass),0.1)*0.7);
    //t_root_depth=fmaxf(2.0, fminf(layer_depth[(nblayers_soil-1)],pow((0.001*total_root_biomass),0.1)*0.7));
    //t_root_depth=fminf(layer_depth[(nblayers_soil-1)], 0.06410256*t_height+1.435897);  // this dependency of tree root depth with tree height is highly questionable -- cf. eg. Stahl et al. 2013 Oecologia, but see also Brum et al. 2019 Journal of Ecology-- but parcimonious. Note that neither Stahl et al. 2013 nor Brum et al. 2019 reported actual root depth, but the "depth of water uptake". Should probably depends more on the total soil depth, and average soil water availabilty -- cf. distribution of maximal root depth per biomes in Canadell/Jackson et al.
    //t_root_depth=0.22*pow(t_dbh*100, 0.54); // this is the allometry used in ED2, Xu et al. 2016, for evergreen trees. Parameter values were fixed based on Kenzo et al. 2009 Journal of Tropical Biology, which reported on data from excavated trees in wet secondary forests in Malaysia, although Xu et al. 2016 applied it in a seasonally dry tropical forests in Costa Rica. This allometry was then compared against the one obtained from excavated trees in a seasonally dry tropical forests in Costa Rica by Smith-Martin et al. 2019 New Phytologist: although the shape and the differences among the allometries for deciduous and evergreen trees were in overall agreement with observations, the ED2 allometry generally underestimated the observed root depth. As it seems that dry forest species have deeper root than wet forest species (see e.g. Holbrook et al. 1995, cited in Smith-Martin et al. 2019, or Markesteijn & Poorter 2009), this should be better suited to our application to wet forest although this remains to be confronted with data. This allometry gives a root depth of 22cm and 2.65 m for trees with dbh=1cm and 1m respectively. Note that in Xu et al., the first parameter (here 0.22) varied among phenological types, so that root depth of evergreen trees is about twice the ones of deciduous. This link between plant phenology and rooting depth was found empirically in Smith-Martin et al. 2019 New Phytologist, Hasselquist et al. 2010 Oecologia, while Markesteijn & Poorter found a correlation between root depth and stem density on first year seedling in Panama. This link between phenology and rooting depth with higher rooting depth for evergreen species was needed to sustain leaf cover in simulations with ED2, despite their higher leaf and stem drought tolerance (P50 and TLP; Smith-Martin et al. 2019). In absence of data, we here assumed this allometry suits to all species but this has to be discussed !
    t_root_depth=fminf(0.35*pow(t_dbh*100, 0.54), layer_depth[(nblayers_soil-1)]); // this is the allometry used in ED2-hydro, Xu et al. 2016, for evergreen trees, with one parameter, b1=root_depth(dbh=1cm), changed to correct the overall underestimation of tree rootdepth it leads to (cf. Smith-Martin et al. 2019 New phytol.). In absence of data on root depth (in general, and even more at the species level, see email with V. Freycon & B. Ferry january 20th 2020), this should probably be calibrated/fine-tuned.
    float i_root_depth=1.0/t_root_depth;
    
    //Tree water potential in the root zone (or predawn water potential)
    float shallow_bound=1.0;
    float deep_bound=1.0;
    float sumG=0.0;
    t_phi_root=0.0;
    
    float layer_depth_previous = 0.0;
    float layer_thickness=0.0;
    float root_area = t_CR*t_CR*PI;

    for(int l=0; l<nblayers_soil; l++) {
        
        layer_thickness=layer_depth[l]- layer_depth_previous;
        layer_depth_previous =layer_depth[l];
        
        deep_bound=exp(-3.0*i_root_depth*layer_depth[l]);
        //deep_bound=exp(-4.02*i_root_depth*layer_depth[l]); // by consistency with the root depth allometry, we here used the parameter value used by Xu et al. 2016 for evergreen trees. This value (log(beta)) was defined based on Jackson et al. 1996 Oecologia.
        //t_soil_layer_weight[l]=shallow_bound-deep_bound;                        //the weight for each soil layer are the relative tree root biomass in this layer. This is the most classic approach, however probably not fully relevant, inasmuch as trees may equilibrate "preferentially" with the wettest part of the soil, where water is more available and more easily exctractable. Hence another possible weighting integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2011). See below.
        t_root_biomass[l]=total_root_biomass*(shallow_bound-deep_bound);            // this is the root biomass in layer l, computed from the integration of the root biomass exponential profile in this layer, following Arora & Boer 2003. This also corresponds to the shape of cumulative root fraction used by Jackson et al. 1996 Oecologia, who gathered a large global database of root distributions with this equation. It was thus used in many models, e.g. Duursma & Medlyn 2012; Xu et al. 2016.
        //t_soil_layer_weight[l]=t_root_biomass[l]*10.0/(-log(sqrt(PI*t_root_biomass[l]*10.0)*0.001)); // this soil layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)
        
        if(t_LA > 0.0){

if (_WATER_RETENTION_CURVE==1) {
            if (_SOIL_LAYER_WEIGHT==0) { // soil layer weights as a function of root biomass in each layer only (cf. M1 in de Kauwe et al 2015)
                t_soil_layer_weight[l]=t_root_biomass[l];
                t_phi_root+=t_soil_layer_weight[l]*soil_phi3D[l][site_DCELL[t_site]];
            } else if (_SOIL_LAYER_WEIGHT==1) {
                t_soil_layer_weight[l]=t_root_biomass[l]*10.0*Ks[l][site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0)*0.001))); // this soil layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)
                t_phi_root+=t_soil_layer_weight[l]*soil_phi3D[l][site_DCELL[t_site]];   //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil layer in the DCELL where the tree stands. Note that KsPhi was here not computed as Ks[l][site_DCELL[t_site]]*soil_phi3D[l][site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                
            } else if (_SOIL_LAYER_WEIGHT==2) { // soil layer weight as in Duursma & Medlyn 2012 (cf M3 in de Kauwe et al. 2015)
                
                //if (soil_phi3D[l][site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=t_root_biomass[l]*10.0*Ks[l][site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0)*0.001)))*(soil_phi3D[l][site_DCELL[t_site]]+3); // this soil layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012). Note that in their implementation of MAESPA, Christina et al. used minimum root water potential = -1.6 MPa and not -3 MPa as here, and they also added a gravimetric component, since they explore the effect of very deep root. Sensibility to this value and addition to be tested.
                //if (soil_phi3D[l][site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=t_root_biomass[l]*10.0*Ks[l][site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0/(length_dcell*length_dcell*layer_thickness))*0.001)))*(soil_phi3D[l][site_DCELL[t_site]]+3);
                if (soil_phi3D[l][site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=(t_root_biomass[l]*10.0/root_area)*Ks[l][site_DCELL[t_site]]/(abs(log(sqrt(PI*t_root_biomass[l]*10.0/(root_area*layer_thickness))*0.001)))*(soil_phi3D[l][site_DCELL[t_site]]+3);
                else t_soil_layer_weight[l]=0.0;
                t_phi_root+=t_soil_layer_weight[l]*soil_phi3D[l][site_DCELL[t_site]];   //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil layer in the DCELL where the tree stands. Note that KsPhi was here not computed as Ks[l][site_DCELL[t_site]]*soil_phi3D[l][site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                
            }
} else if (_WATER_RETENTION_CURVE==0) {
    
            if (_SOIL_LAYER_WEIGHT==0) { // soil layer weights as a function of root biomass in each layer only (cf. M1 in de Kauwe et al 2015)
                t_soil_layer_weight[l]=t_root_biomass[l];
                t_phi_root+=t_soil_layer_weight[l]*soil_phi3D[l][site_DCELL[t_site]];
            } else if (_SOIL_LAYER_WEIGHT==1) {
                t_soil_layer_weight[l]=t_root_biomass[l]*10.0/(abs(log(sqrt(PI*t_root_biomass[l]*10.0)*0.001))); // this soil layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012)
                t_phi_root+=t_soil_layer_weight[l]*KsPhi[l][site_DCELL[t_site]];   //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil layer in the DCELL where the tree stands. Note that KsPhi was here not computed as Ks[l][site_DCELL[t_site]]*soil_phi3D[l][site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                t_soil_layer_weight[l]*=Ks[l][site_DCELL[t_site]];
                
            } else if (_SOIL_LAYER_WEIGHT==2) { // soil layer weight as in Duursma & Medlyn 2012 (cf M3 in de Kauwe et al. 2015)

                if (soil_phi3D[l][site_DCELL[t_site]]>(-3.0)) t_soil_layer_weight[l]=(t_root_biomass[l]*10.0/root_area)/(abs(log(sqrt(PI*t_root_biomass[l]*10.0/(root_area*layer_thickness))*0.001)))*(soil_phi3D[l][site_DCELL[t_site]]+3); // this soil layer weight integrates the soil-to-root conductance into account (as in de Kauwe et al. 2015; Duursma & Medlyn 2012). Note that in their implementation of MAESPA, Christina et al. used minimum root water potential = -1.6 MPa and not -3 MPa as here, and they also added a gravimetric component, since they explore the effect of very deep root. Sensibility to this value and addition to be tested.
                else t_soil_layer_weight[l]=0.0;
                t_phi_root+=t_soil_layer_weight[l]*KsPhi[l][site_DCELL[t_site]];   //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil layer in the DCELL where the tree stands. Note that KsPhi was here not computed as Ks[l][site_DCELL[t_site]]*soil_phi3D[l][site_DCELL[t_site]], to avoid some potential divergence for very low water content, and due to the limit of float type, but directly as the exact power of SWC -- see in UpdateField.
                t_soil_layer_weight[l]*=Ks[l][site_DCELL[t_site]];
                
            }
            
}
            
            
            if(isnan(t_phi_root)) {
                cout << "nan t_phi_root: " << endl;
                cout << l << "\t" << t_soil_layer_weight[l] << "\t" << KsPhi[l][site_DCELL[t_site]] << soil_phi3D[l][site_DCELL[t_site]] << "\t" << endl;
            }
            
        } else {
            if (l==0) {
                t_soil_layer_weight[l]=1;
                t_phi_root+=soil_phi3D[l][site_DCELL[t_site]];
            }
            else t_soil_layer_weight[l]=0;
        }
        
        //if(t_soil_layer_weight[l]<=0) cout << "t_soil_layer_weight[l]=" << t_soil_layer_weight[l] << " t_root_biomass[l]=" << t_root_biomass[l] << " Ks=" << Ks[l][site_DCELL[t_site]] << " Ks*Phi soil =" << KsPhi[l][site_DCELL[t_site]] << " phi soil=" <<   soil_phi3D[l][site_DCELL[t_site]] << endl;
        //t_soil_layer_weight[l]=t_root_biomass[l]*10/(-log(sqrt(PI*t_root_biomass[l]*10)*0.001));
        //t_phi_root+=t_soil_layer_weight[l]*fmaxf(0.0,(KsPhi2[l][site_DCELL[t_site]]-KsPhi[l][site_DCELL[t_site]]*t_s->s_tlp));  //the water potential in the root zone is computed as the weighted mean of the soil water potential in each soil layer
        //t_soil_layer_weight[l]*=fmaxf(0.0,(KsPhi[l][site_DCELL[t_site]]-Ks[l][site_DCELL[t_site]]*t_s->s_tlp));
        sumG+=t_soil_layer_weight[l];
        shallow_bound=deep_bound;
        
        if(t_soil_layer_weight[l]<0.0 || isnan(t_soil_layer_weight[l])) {
            cout << "incorrect soil_layer_weight: " << endl;
            cout << l <<  "\t" << site_DCELL[t_site]<<  "\t" << t_soil_layer_weight[l] << "\t" << t_root_biomass[l] << "\t" << Ks[l][site_DCELL[t_site]] << "\t" << KsPhi[l][site_DCELL[t_site]] << "\t" << -log(sqrt(PI*t_root_biomass[l]*10)*0.001) << "\t" << deep_bound << "\t" << shallow_bound << "\t" << SWC3D[l][site_DCELL[t_site]] <<  "\t" <<(SWC3D[l][site_DCELL[t_site]]-Min_SWC[l])/(Max_SWC[l]-Min_SWC[l]) <<  "\t"  << Ksat[l]  <<  "\t" << soil_phi3D[l][site_DCELL[t_site]];
            cout <<endl;
        }
        
    }
    
    
    if (sumG>0.0) {
        float isumG=1.0/sumG;
        for (int l=0; l<nblayers_soil; l++) {
            t_soil_layer_weight[l]*=isumG;
        }
        t_phi_root*=isumG;
        if (isnan(t_phi_root)) {
            cout << "nan t_phi_root" << endl;
            cout << isumG << endl;
        }
    } else {
        float iRootB=1.0/total_root_biomass;
        for (int l=0; l<nblayers_soil; l++) {
            t_phi_root+=t_root_biomass[l]*soil_phi3D[l][site_DCELL[t_site]];
            if (isnan(t_phi_root)) {
                cout << "nan t_phi_root" << endl;
                cout << t_root_biomass[l] << "\t" << soil_phi3D[l][site_DCELL[t_site]] << endl;
            }
            t_soil_layer_weight[l]=t_root_biomass[l]*iRootB;
        }
        t_phi_root*=iRootB;
        if (isnan(t_phi_root)) {
            cout << "nan t_phi_root" << endl;
            cout << iRootB << endl;
        }
    }
    
    float sum_weight=0.0;
    for (int l=0; l<nblayers_soil; l++) {
        sum_weight+=t_soil_layer_weight[l];
    }
    if (sum_weight>1.5) {
        cout << "Weird weights in Water_availability; sum_weight=" << sum_weight << " sumG=" << sumG << " total_root_biomass=" << total_root_biomass << " t_phi_root=" << t_phi_root ;
        if (sumG>0.0) {
            float isumG=1.0/sumG;
            cout << " isumG=" << isumG << endl;
        }
        cout << "t_soil_layer_weight[l]: ";
        for (int l=0; l<nblayers_soil; l++) {
            cout <<t_soil_layer_weight[l] << "\t";
        }
        cout << endl;
        cout << "soil_phi: ";
        for (int l=0; l<nblayers_soil; l++) {
            cout << soil_phi3D[l][site_DCELL[t_site]] << "\t";
        }
        cout << endl;
        
    }
    
    t_phi_root-=0.01*t_height; // gravitational effect on leaf water potential and stress
    
    if (isnan(t_phi_root)) {
        cout << "nan t_phi_root" << endl;
        cout << t_height << endl;
    }
    
    //! Tree water stress factors
    
    //t_WSF=fminf(1.0, fmaxf(0.0, ((SWC3D[0][site_DCELL[t_site]]-Min_SWC[0])/(Max_SWC[0]-Min_SWC[0])))); // this is the simple linear WSF, with SWC as independent varible, used in lots of model (see Powell et al. 2013 New Phytol, De Kauwe et al. 2015 Biogeosciences, Laio et al. 2001 Advances in Water resources, Egea et al. 2011 AFM....)
    //t_WSF=fminf(1.0, fmaxf(0.0, (t_phi_root-(t_s->s_tlp))*(t_s->s_dWSF))); // this is the linear WSF, with wtare potential as independent variable, used in CLM model (see Powell et al. 2013 New Phytologist, Verhoef & Egea 2011 AFM)
    //t_WSF=exp((log(0.05)/t_s->s_tlp)*t_phi_root); // this is the WSF, to simulate stomatal limitation (ie. hinder g1), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences. If this version of WSF is finally adopted, the parameter b=log(0.05)/t_s->s_tlp, should be declare as species variable (instead of s_dWSF).
    t_WSF=exp(t_b*t_phi_root); // this is the WSF shape used to simulate stomatal limitation (ie. hinder g1), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences, but with a different parameterization for the parameter b (using the relationship between phi_gs90 and tlp from Martin-StPaul et al. 2017 Ecology letters, and assuming the WSF=0.9 at phi_gs90).
    //float par=t_s->s_tlp+1;
    //t_WSF_A=(1.0+exp(6.0*par))/(1.0+exp(6.0*(par-t_phi_root))); // this is the WSF, to simulate non-stomatal limitation (ie. hinder Vcmax and Jmax), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences
    //t_WSF_A=(1.0+exp(3.0*par))/(1.0+exp(3.0*(par-t_phi_root))); // this is the WSF, to simulate non-stomatal limitation (ie. hinder Vcmax and Jmax), drawn from Zhou et al. 2013 AFM, and de Kauwe et al. 2015 Biogeosciences
#ifdef WATER
    //t_g1= (-3.97 * t_wsg + 6.53)*t_WSF; this is the relationship provided by Lin et al. 2015
    t_g1 = t_g1_0*t_WSF; // with the water stress factor added
#endif
    
    t_WSF_A=1/(1+pow(t_phi_root*t_itlp, 6)); // this is the WSF used to hinder Vcmax and Jmax in Xu et al. 2016 (equ. S5) (ED2-hydro).
    
    
    if (t_WSF < 0.0 || t_WSF_A < 0.0 || t_WSF>1.0 || t_WSF_A >1.0 ||t_phi_root >0.0 || isnan(t_WSF) || isnan(t_WSF_A) || isnan(t_phi_root)) {
        cout << "incorrect value in one of WSF, WSF_A, t_phi_root " << endl;
        cout <<t_WSF << "\t" << t_phi_root << "\t" << t_tlp  << "\t" << t_dbh << "\t" << t_height << "\t" << t_age << "\t" << t_phi_lethal << "\t" << sumG;
        cout <<endl;
    }
    //}
}

//###################################################################
// Contribution of trees to the stand Transpiration field. Called by UpdateField
//####################################################################
//! - Adds up each tree contribution to the stand Transpiration field, that is the water removed from the soil through all tree transpiration.
//! - Similar to CalcLAI, that adds up each tree contribution to LAI3D field.
void Tree::Water_uptake() {
    if(t_age > 0.0) {
        int l = 0;
        float depth = 0.0;
        float sum_weights=0.0;
        while (sum_weights<1.0 && l < nblayers_soil) {
            Transpiration[l][site_DCELL[t_site]] += t_soil_layer_weight[l] * t_transpiration;
            if (t_soil_layer_weight[l] < 0.0 || t_transpiration < 0.0 ) {
                cout << setprecision(10);
                cout << "Problem with soil_layer_weight and transpiration at site: " << t_site << " layer: " << l << " depth: " << depth << endl;
                cout << t_soil_layer_weight[l] << "\t" << t_transpiration << "\n";
            }
            sum_weights+=t_soil_layer_weight[l];
            depth = layer_depth[l];
            l++;
        }
    }
}

#endif

                // update derived parameters
                sites = rows*cols;
                sites_per_dcell = length_dcell*length_dcell;
                nbdcells = int(sites/sites_per_dcell);
                linear_nb_dcells = int(cols/length_dcell);
                cout << "rows: " << rows << " cols: " << cols << " HEIGHT: " << HEIGHT << endl;
                cout << "Number of dcells: " << nbdcells << endl;
                cout << "Lin number of dcells: " << linear_nb_dcells << endl;
                cout << _WATER_RETENTION_CURVE << " " << _SOIL_LAYER_WEIGHT << " " << Cair << endl;
                
#ifdef WATER
                i_sites_per_dcell=1.0/float(sites_per_dcell);
                PPFDtoSW= 1/SWtoPPFD;

#ifdef WATER
        //! Global function: This function reads inputs from the soil input file
        //! - in v.3.0, all soil parameters (Sat_SWC, Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (Sat_SWC, Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (Sat_SWC, Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
        //! - Sat_SWC and Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
        void ReadInputSoil(){
            cout << endl << "Reading in file: " << inputfile_soil << endl;
            fstream InSoil(inputfile_soil, ios::in);
            
            if(InSoil){
                InSoil.getline(buffer,256,'\n');
                vector<float> layer_thickness, proportion_Silt, proportion_Clay, proportion_Sand; // in m, %, %,%
                vector<float> SOC, DBD, pH, CEC; //soil organic content, provided in %; dry bulk density, in g cm-3; pH; cation exchange capacity, in cmol kg-1
                SOC.reserve(20);
                DBD.reserve(20);
                pH.reserve(20);
                CEC.reserve(20);

                layer_thickness.reserve(20);
                proportion_Silt.reserve(20);
                proportion_Clay.reserve(20);
                proportion_Sand.reserve(20);
                nblayers_soil = 0;
                
                // we go through all lines in the input file
                string line;
                while(getline(InSoil, line)){
                    istringstream linestream(line);
                    
if (_WATER_RETENTION_CURVE==1) {
                    float thickness_current, proportion_Silt_current, proportion_Clay_current, proportion_Sand_current, SOC_current, DBD_current, pH_current, CEC_current;
                    linestream >> thickness_current >> proportion_Silt_current >> proportion_Clay_current >> proportion_Sand_current >> SOC_current >> DBD_current >> pH_current >> CEC_current;
                    layer_thickness.push_back(thickness_current);
                    proportion_Silt.push_back(proportion_Silt_current);
                    proportion_Clay.push_back(proportion_Clay_current);
                    proportion_Sand.push_back(proportion_Sand_current);
                    SOC_current*=10; // to convert SOC in % to gC kg-1 soil, to apply Hodnett & Tomasella 2002
                    SOC.push_back(SOC_current);
                    DBD.push_back(DBD_current);
                    pH.push_back(pH_current);
                    CEC.push_back(CEC_current);
} else if (_WATER_RETENTION_CURVE==0) {
                    float thickness_current, proportion_Silt_current, proportion_Clay_current, proportion_Sand_current;
                    linestream >> thickness_current >> proportion_Silt_current >> proportion_Clay_current >> proportion_Sand_current;
                    layer_thickness.push_back(thickness_current);
                    proportion_Silt.push_back(proportion_Silt_current);
                    proportion_Clay.push_back(proportion_Clay_current);
                    proportion_Sand.push_back(proportion_Sand_current);
}
                    
                    nblayers_soil++;
                }
                cout << "Read in: " << nblayers_soil << " soil layers" << endl;
                
                if(NULL==(layer_depth=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc layer_depth" << endl;
                float cumulative_depth=0.0;
                for (int l=0; l<nblayers_soil; l++) {
                    cumulative_depth+=layer_thickness[l];
                    layer_depth[l]=cumulative_depth;
                }
                // (added in the header but kept in here) in this version, all soil parameters (Sat_SWC, Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (Sat_SWC, Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (Sat_SWC, Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
                // (added in the header but kept in here) Sat_SWC and Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
                if(NULL==(Sat_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Sat_SW" << endl;
                if(NULL==(Max_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Max_SW" << endl;
                if(NULL==(FC_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc FC_SW" << endl;

                
                for (int l=0; l<nblayers_soil; l++) {
                    
if (_WATER_RETENTION_CURVE==1) {
                    Sat_SWC[l] = 0.01*(81.799+(0.099*proportion_Clay[l])-(31.42*DBD[l])+(0.018*CEC[l])+(0.451*pH[l])-(0.0005*proportion_Sand[l]*proportion_Clay[l]));  // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3
} else if (_WATER_RETENTION_CURVE==0) {
                    Sat_SWC[l]= 0.01*(40.61+(0.165*proportion_Silt[l])+(0.162*proportion_Clay[l])+(0.00137*proportion_Silt[l]*proportion_Silt[l])+(0.000018*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3
}
                    Max_SWC[l]=Sat_SWC[l]*sites_per_dcell*LH*LH*layer_thickness[l]; // in m3
                    cout << "layer " << l << " Vol=" << sites_per_dcell*LH*LH*layer_thickness[l]<< " m3; Sat_SWC =" << Sat_SWC[l] << " MAX_SWC =" << Max_SWC[l] << " m3." << endl;
                }
                
                if(NULL==(Ksat=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Ksat" << endl;
                for (int l=0; l<nblayers_soil; l++) {
                    Ksat[l]=0.007055556*pow(10,(-0.60-(0.0064*proportion_Clay[l])+(0.0126*proportion_Sand[l]))); // according to Cosby et al. 1984 (the only expression of k_sat reported in Table 2 of Marthews et al. 2014). k_sat is here in mm/s or equivalently in kg/m2/s.
                    cout << "layer " << l << " Ksat=" << Ksat[l] <<"mm/s or kg/m2/s  "<< Ksat[l]*9.8/18 << endl;
                }
                
                if(NULL==(Res_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Res_SW" << endl;
                if(NULL==(Min_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Max_SW" << endl;
                
                for (int l=0; l<nblayers_soil; l++) {
                    
if (_WATER_RETENTION_CURVE==1) {
                    Res_SWC[l]= 0.01*(22.733-(0.164*proportion_Sand[l])+(0.235*CEC[l])-(0.831*pH[l])+(0.0018*proportion_Clay[l]*proportion_Clay[l])+(0.0026*proportion_Sand[l]*proportion_Clay[l])); // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3
} if (_WATER_RETENTION_CURVE==0) {
                    Res_SWC[l]= 0.01*fmaxf(0.0,(-2.094+(0.047*proportion_Silt[l])+(0.431*proportion_Clay[l])-(0.00827*proportion_Silt[l]*proportion_Clay[l]))); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
}
                    Min_SWC[l]=Res_SWC[l]*sites_per_dcell*LH*LH*layer_thickness[l]; //in m3

                    cout << "layer " << l << " Vol=" << sites_per_dcell*LH*LH*layer_thickness[l] << "m3; Res=" << Res_SWC[l]<<  " MIN_SWC =" << Min_SWC[l] << " m3" << endl;
                }
                
if (_WATER_RETENTION_CURVE==1) {
                if(NULL==(a_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc a_vgm" << endl;
                if(NULL==(b_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc b_vgm" << endl;
                if(NULL==(c_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc c_vgm" << endl;
                if(NULL==(m_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc m_vgm" << endl;

                for (int l=0; l<nblayers_soil; l++) {

                    float alpha=1000*exp((-2.294-(3.526*proportion_Silt[l])+(2.440*(0.1*SOC[l]))-(0.076*CEC[l])-(11.331*pH[l])+(0.019*proportion_Silt[l]*proportion_Silt[l]))*0.01); // this is the alpha parameter of the van Genuchten-Mualem model, in MPa-1 (ie after already dividing by rho*g, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                    a_vgm[l]=-1.0/alpha;
                    float n_vgm=exp((62.986-(0.833*proportion_Clay[l])-(0.529*(SOC[l]*0.1))+(0.593*pH[l])+(0.007*proportion_Clay[l]*proportion_Clay[l])-(0.014*proportion_Sand[l]*proportion_Silt[l]))*0.01);    // this is the n parameter of the van Genuchten-Mualem model, dimensionless, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                    m_vgm[l]=1.0-1.0/n_vgm;
                    b_vgm[l]=1.0/m_vgm[l];
                    c_vgm[l]=1.0-m_vgm[l];
                    
                    FC_SWC[l]=(Res_SWC[l] + (Sat_SWC[l]-Res_SWC[l])*pow((pow(0.01*alpha, 1/c_vgm[l])+1),-(1/b_vgm[l])))*sites_per_dcell*LH*LH*layer_thickness[l]; // this is the layer water content at field capacity, in m3. As in Marthews et al. 2014 (cf. note in Table 2), we used Phi at FC=-10kPa and not -33kPa, following Marshall et al., 1996; Townend et al., 2001; Tomasella and Hodnett, 2004)
                    
                    cout << "layer " << l << " alpha=" << alpha << "\t" << "n_vgm=" << n_vgm << " FC_SWC=" << FC_SWC[l] << endl;
                }
} else if (_WATER_RETENTION_CURVE==0) {
                if(NULL==(phi_e=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc phi_e" << endl;
                if(NULL==(b=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc b" << endl;
                for (int l=0; l<nblayers_soil; l++) {
                    phi_e[l]=-0.001*(0.285+(0.000733*proportion_Silt[l]*proportion_Silt[l])-(0.00013*proportion_Silt[l]*proportion_Clay[l])+(0.0000036*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. phi_e is here in MPa.
                    b[l]=exp(1.197+(0.00417*proportion_Silt[l])-(0.0045*proportion_Clay[l])+(0.000894*proportion_Silt[l]*proportion_Clay[l])-(0.00001*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
                    //phi_e[l]=-0.00000001*pow(10.0,(2.17-(0.0063*proportion_Clay[l])-(0.0158*proportion_Sand[l])))*(1000*9.80665); // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. In MPa.
                    //b[l]=3.10+0.157*proportion_Clay[l]-0.003*proportion_Sand[l]; // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. Dimensionless.
                    
                    FC_SWC[l]=(Res_SWC[l] + (Sat_SWC[l]-Res_SWC[l])*pow(-0.01/phi_e[l],-(1/b[l])))*sites_per_dcell*LH*LH*layer_thickness[l]; // this is the layer water content at filed capacity, in m3. As in Marthews et al. 2014 (cf. note in Table 2), we used Phi at FC=-10kPa and not -33kPa, following Marshall et al., 1996; Townend et al., 2001; Tomasella and Hodnett, 2004)
                    
                    
                    cout << "layer " << l << " phi_e=" << phi_e[l] << "\t" << "b=" << b[l] << " FC_SWC=" << FC_SWC[l] << endl;
                }
}
                cout << "Successfully read the soil file" << endl;
                
                
            } else {
                cout<< "ERROR with the soil file" << endl;
            }
            cout << endl;
        }
#endif
        
        //! Global function: initialise intraspecific variables
        void InitialiseIntraspecific(){
            float max_intraspecific_height=0.0, min_intraspecific_height=1000.0,
            max_intraspecific_CR=0.0, min_intraspecific_CR = 1000.0,
            max_intraspecific_CD=0.0, min_intraspecific_CD = 1000.0,
            max_intraspecific_P=0.0, min_intraspecific_P = 1000.0,
            max_intraspecific_N=0.0, min_intraspecific_N = 1000.0,
            max_intraspecific_LMA=0.0, min_intraspecific_LMA = 1000.0,
            max_intraspecific_wsg=0.0, min_intraspecific_wsg = 1000.0,
            max_intraspecific_dbhmax=0.0, min_intraspecific_dbhmax = 1000.0;
#ifdef WATER
            float max_intraspecific_leafarea=0.0, min_intraspecific_leafarea = 1000.0,
            max_intraspecific_tlp=0.0, min_intraspecific_tlp = 1000.0;
            double variation_leafarea=0.0,variation_tlp=0.0;
#endif
            
            double variation_height=0.0, variation_CR=0.0, variation_CD=0.0, variation_P=0.0, variation_N=0.0, variation_LMA=0.0,variation_wsg=0.0,variation_dbhmax=0.0;
            for(int i=0;i<10000;i++){ // modified FF v.3.1.5 (reduced from 100000 to 10000)
                if(covariance_status == 0){
                    variation_N = gsl_ran_gaussian(gslrand, sigma_N);
                    variation_P = gsl_ran_gaussian(gslrand, sigma_P);
                    variation_LMA = gsl_ran_gaussian(gslrand, sigma_LMA);
                } else {
                    gsl_ran_multivariate_gaussian(gslrand, mu_N_P_LMA, mcov_N_P_LMA, variation_N_P_LMA);
                    variation_N = gsl_vector_get(variation_N_P_LMA, 0);
                    variation_P = gsl_vector_get(variation_N_P_LMA, 1);
                    variation_LMA = gsl_vector_get(variation_N_P_LMA, 2);
                }
                gsl_ran_bivariate_gaussian(gslrand, sigma_height, sigma_CR, corr_CR_height, &variation_height, &variation_CR);
                // variation_height = gsl_ran_gaussian(gslrand, sigma_height);
                // variation_CR = gsl_ran_gaussian(gslrand, sigma_CR);
                variation_CD = gsl_ran_gaussian(gslrand, sigma_CD);
                variation_wsg = gsl_ran_gaussian(gslrand, sigma_wsg);
                variation_dbhmax = gsl_ran_gaussian(gslrand, sigma_dbhmax);
#ifdef WATER
                variation_leafarea = gsl_ran_gaussian(gslrand, sigma_leafarea);
                variation_tlp = gsl_ran_gaussian(gslrand, sigma_tlp);
#endif

if (_WATER_RETENTION_CURVE==1) {
                    output[33] << "layer" << "\t" << "depth" << "\t" << "sat" << "\t" << "max" << "\t" << "fc" << "\t" << "res" << "\t" << "min" << "\t" << "Ksat" << "\t" << "a_vgm" << "\t" << "m_vgm" << endl;
                    for (int l=0; l<nblayers_soil; l++) {
                    output[33] << l  << "\t" << layer_depth[l]  << "\t" << Sat_SWC[l] << "\t" << Max_SWC[l]  << "\t" << FC_SWC[l] << "\t" << Res_SWC[l] << "\t" << Min_SWC[l] << "\t" << Ksat[l] << "\t" << a_vgm[l] << "\t" << m_vgm[l] << endl;
                    }

} else if (_WATER_RETENTION_CURVE==0) {
                    
                    output[33] << "layer" << "\t" << "depth" << "\t" << "sat" << "\t" << "max" << "\t" << "fc" << "\t" << "res" << "\t" << "min" << "\t" << "Ksat" << "\t" << "phi_e" << "\t" << "b" << endl;
                    for (int l=0; l<nblayers_soil; l++) {
                    output[33] << l  << "\t" << layer_depth[l]  << "\t" << Sat_SWC[l] << "\t" << Max_SWC[l] << "\t" << FC_SWC[l] << "\t" << Res_SWC[l] << "\t" << Min_SWC[l] << "\t" << Ksat[l] << "\t" << phi_e[l] << "\t" << b[l] << endl;
                         }
}
                   
#endif

#ifdef WATER
                
                for(int l=0;l<nblayers_soil;l++){
                    char par_name[30];
                    sprintf(par_name, "root_biomass%d", l);
                    parameter_names.push_back(par_name);
                }
                nb_parameters += nblayers_soil;
                for(int l=0;l<nblayers_soil;l++){
                    char par_name[30];
                    sprintf(par_name, "soil_layer_weight%d", l);
                    parameter_names.push_back(par_name);
                }
                nb_parameters += nblayers_soil;
                
#endif

            // compute soil_phi3D for RecruitTree function
            for (int d=0; d<nbdcells; d++) {
                for (int l=0; l<nblayers_soil; l++) {
                    float theta_w=(SWC3D[l][d]-Min_SWC[l])/(Max_SWC[l]-Min_SWC[l]);
                if (_WATER_RETENTION_CURVE==1) {
                    if(theta_w==0) {
                        theta_w=0.001; // SS addition for limit value
                        cout << "Warning theta_w = 0 " << endl ;
                    }
                    soil_phi3D[l][d]=a_vgm[l]*pow((pow(theta_w,-b_vgm[l])-1), c_vgm[l]); // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                    float inter= 1-pow((1-pow(theta_w, b_vgm[l])),m_vgm[l]);
                    Ks[l][d]=Ksat[l]*pow(theta_w, 0.5)*inter*inter; // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                    if (isnan(soil_phi3D[l][d]) || isnan(Ks[l][d]) ||  (SWC3D[l][d]-Min_SWC[l])<0) //|| KsPhi[l][d]==0.0 || Ks[l][d]==0.0 || soil_phi3D[l][d]==0.0)
                        cout << "In bucket model, layer " << l << " dcell " << d << " theta_w=" << theta_w << " SWC3D[l][d]-Min_SWC[l]=" << (SWC3D[l][d]-Min_SWC[l]) << " soil_phi3D[l][d]=" << soil_phi3D[l][d] << " Ksat=" << Ksat[l] << " Ks[l][d]=" << Ks[l][d] << endl ;
                } else if (_WATER_RETENTION_CURVE==0) {
                    soil_phi3D[l][d]=phi_e[l]*pow(theta_w, -b[l]); // this is the soil water characteristic of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                    Ks[l][d]=Ksat[l]*pow(theta_w, 2.5+2*b[l]); // this is the hydraulic conductivity curve of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                    KsPhi[l][d]=Ksat[l]*phi_e[l]*pow(theta_w, 2.5+b[l]); //Ks times soil_phi3D, computed directly as the exact power of theta.
                    if (isnan(soil_phi3D[l][d]) || isnan(Ks[l][d]) ||  isnan(KsPhi[l][d]) || (SWC3D[l][d]-Min_SWC[l])<0) //|| KsPhi[l][d]==0.0 || Ks[l][d]==0.0 || soil_phi3D[l][d]==0.0)
                        cout << "In bucket model, layer " << l << " dcell " << d << " theta_w=" << theta_w << " SWC3D[l][d]-Min_SWC[l]=" << (SWC3D[l][d]-Min_SWC[l]) << " soil_phi3D[l][d]=" << soil_phi3D[l][d] << " Ksat=" << Ksat[l] << " phi_e=" << phi_e[l] <<" b[l]=" << b[l] << " KsPhi[l][d]=" << KsPhi[l][d] << " Ks[l][d]=" << Ks[l][d] << endl ;
                 }
                }
            }

            // compute LAID for RecruitTree function
            for(int h=0;h<(HEIGHT+1);h++)
                for(int sbsite=0;sbsite<sites+2*SBORD;sbsite++)
                    LAI3D[h][sbsite] = 0.0;
            for(int site=0;site<sites;site++) T[site].CalcLAI(); // Each tree contribues to LAI3D
            for(int h=HEIGHT;h>0;h--){ // LAI is computed by summing LAI from the canopy top to the ground
                for(int site=0;site<sites;site++){
                    int sbsite=site+SBORD;
                    LAI3D[h-1][sbsite] += LAI3D[h][sbsite];
                }
            }
            
#ifdef WATER
            for (int d=0; d<nbdcells; d++) {
                Runoff[d]=0.0;
                Interception[d]=0.0;
                Throughfall[d]=0.0;
                Evaporation[d]=0.0;
                Leakage[d]=0.0;
                //for (int l=0;l<nblayers_soil;l++) {
                //    Transpiration[l][d]=0.0;
                //}
                for(int h=0;h<(HEIGHT+1);h++) {
                    LAI_DCELL[h][d]=0.0;
                }
                Canopy_height_DCELL[d]=0.0;
                HSum_DCELL[d]=0;
                TopWindSpeed_DCELL[d]=0.0;
            }

#ifdef WATER
            //**  Evolution of belowground hydraulic fields: Soil bucket model
            
            //for(int site=0;site<sites;site++) T[site].Water_uptake(); // Update of Transpiration: tree water uptake, each tree will deplete soil water content through its transpiration. Now made ate the end of the evolution loop so that the outputs for water uptake match the others (otherwise lag of one timestep)
            
            for (int d=0; d<nbdcells; d++) {
                //****   BUCKET MODEL in each dcell   ****
                // the unit used for water volume throughout the bucket model is m3.
                // NOTE: under the assumption of a flat terrain and no lateral fluxes, as it is assumed here for a first implementation, the order with which dcells are visited during the loop does not matter. With topography, we will need to visit the soil voxels (ie. dcells*layer) from highest to lowest elevation so that run-off from highest voxels contribute to the water flux entering the lowest ones.
                //  to be investigated: does the order in which transpiration and evaporation are retrieved from the soil affect the overall outcome? which one should be retrieved first?
                
                

                
                // Water uptake through tree transpiration
                float w_uptake=0.0;
                for (int l=0; l<nblayers_soil; l++) {
                    w_uptake=fminf(Transpiration[l][d],(SWC3D[l][d]-Min_SWC[l]));
                    if (Transpiration[l][d]<0.0 ||  isnan(Transpiration[l][d]) || isnan(w_uptake) || (SWC3D[l][d]-Min_SWC[l])<0 ) {
                        cout << "l=" << l << " d= "<< d <<" transpiration=" << Transpiration[l][d] << " and SWC3D[l][d]=" << SWC3D[l][d] << " and Min_SWC[l]=" << Min_SWC[l] << " and SWC3D[l][d]-Min_SWC[l]=" << SWC3D[l][d]-Min_SWC[l] << endl;
                    }
                    SWC3D[l][d]-=w_uptake;
                    SWC3D[l][d]=fmaxf(SWC3D[l][d],Min_SWC[l]);
                    if (Transpiration[l][d]<0.0 ||  isnan(Transpiration[l][d]) || isnan(w_uptake) || (SWC3D[l][d]-Min_SWC[l])<0 ) {
                        cout << "After w_uptake, l=" << l << " d= "<< d <<" transpiration=" << Transpiration[l][d] << " and SWC3D[l][d]=" << SWC3D[l][d] << " and Min_SWC[l]=" << Min_SWC[l] << " and SWC3D[l][d]-Min_SWC[l]=" << SWC3D[l][d]-Min_SWC[l] << endl;
                    }
                }
                
                // Evaporation from soil
                
                // if this should be negligible in dense forest understory, it should have a more important effect in open areas, especially through species filtering at germination stage, at the beginning of a succession or in gaps in drier conditions; see Marthews et al. 2008 Ecological Modelling
                // However it is sometimes neglected and not represented in models, eg. Laio et al. 2001, Guterriez et al. 2014, Fischer et al. 2014.
                // note that, in this version, evaporation only depletes the most superficial soil layer. This could be changed, especially if the superficial soil layer is particularly thin and the energy reaching the soil high.
                
                
                // here, we use a phenomenological approach, following Granier et al. 1999 Ecological Modelling and Wagner et al. 2011 AFM, which assumed that evaporation is proportional to the energy reaching the soil.[this is an approximation as as the soil gets drier, more energy would be needed to remove the same amount of water from the soil as water molecules should be more tighly bound to soil particules and cavitation also occur in the soil...] ==> see if a model under which evaporation also depends on the soil water potential would not be better -- I guess so.
                // parameter values are not so clear, so TO BE CHECKED.
                // float e_factor=PPFDtoSW * 3600*0.000001*nbhours_covered* 0.1 * sites_per_dcell*LH*LH*0.001; // to be moved outside of the loop to avoid repeating calculation.
                //float e_Granier = e_factor* WDailyMean * exp(-klight*LAI_DCELL[0][d]);
                //3600*0.000001*nbhours_covered to convert Wmax in micromol of PAR /s /m2 into  Joule, and 10^-6 to MJoule as in Wagner et al. 2011 (however the value provided by Wagner et al. 2011 seems really weird -too high-, and the values we obtained here are in agreement with the ones reported in Marthews et al. 2014.
                //the value 0.1 is drawn from Wagner et al. 2011, but not really explained... to be checked!
                //sites_per_dcell*LH*LH*0.001 is to convert the amount of water in mm, ie. in 10-3 m3/m2, to the amount of water evaporated for the focal dcell in m3
                
                // in this newer version, we used the framework provided by Sellers et al. 1992, which is better mechanistically grounded: depends on the soil layer resistance, which varies with its water potential, and the aerodynamic resistance in series and the differences of vapour pressure between the top soil layer and air just above
                float absorb_prev = LAI_DCELL[1][d];
                float absorb_current = LAI3D[0][d];
                float absorb_delta = absorb_current - absorb_prev;
                if(absorb_delta < 0.0) absorb_delta = 0.0;    //! eliminate rounding errors
                int intabsorb = CalcIntabsorb(absorb_prev, absorb_delta);
                float VPDground = VPDDailyMean*LookUp_VPD[intabsorb]*1000; // in Pa
                float Tsoil = tDailyMean - LookUp_T[intabsorb];
                float esat_ground= 611.21*exp((18.678 - (Tsoil/234.5))*(Tsoil/(257.14+Tsoil))); // Buck equation; in Pa (see Jones p. 348)
                float esoil= esat_ground*exp(2.17*soil_phi3D[0][d]/(Tsoil-ABSZERO)); // esoil variation with the top soil layer water potential, following Duursma & Medlyn 2012 equ. 17, Cochard et al. 2021 equ. 36., see equ. 5.14 in Jones (p. 102), in Pa
                float eair= esat_ground-VPDground; // in Pa
                //float r_soil = exp(8.206 - 4.255*SWC3D[0][d]/Max_SWC[0]) ; // soil surface resistance in s m-1, following Sellers et al. 1992 equ. 19, see also equ 12 in Merlin et al. 2016 (also used in CLM, Oleson et al. 2007).
                float r_soil = exp(8.206 - 4.255*SWC3D[0][d]/FC_SWC[0]) ; // soil surface resistance in s m-1, following Sellers et al. 1992 equ. 19, see also equ 12 in Merlin et al. 2016 (also used in CLM, Oleson et al. 2007).
#ifdef FULL_CLIMATE
                float r_aero = 43.17347*exp(alphaInoue*(1-1/Canopy_height_DCELL[d]))/(WSDailyMean*TopWindSpeed_DCELL[d]); // aerodynamic resistance to hear transfer (boundary layer just above the soil surface), in s m-1 (see equ. 7 and 14 in Duursma & Medlyn 2012; and equ. B10 in Merlin et al. 2016). 43.17347= log(1/0.001)/(0.40*0.40), where 1= the reference height where the wind speed is measured, in m, 0.001=the momentum soil roughness in m (set to 0.001 following Yang et al. 2008 and Stefan et al 2015 in Merlin et al. 2016 equ B10), and 0.40=the von Karman constant.
#else
                float r_aero = 43.17347*exp(alphaInoue*(1-1/Canopy_height_DCELL[d]))/TopWindSpeed_DCELL[d]; // aerodynamic resistance to hear transfer (boundary layer just above the soil surface), in s m-1 (see equ. 7 and 14 in Duursma & Medlyn 2012; and equ. B10 in Merlin et al. 2016). 43.17347= log(1/0.001)/(0.40*0.40), where 1= the reference height where the wind speed is measured, in m, 0.001=the momentum soil roughness in m (set to 0.001 following Yang et al. 2008 and Stefan et al 2015 in Merlin et al. 2016 equ B10), and 0.40=the von Karman constant.
#endif
                float Rtot= r_soil+r_aero; // in s m-1
                float e = nbhours_covered*sites_per_dcell*LH*LH*0.0078*(esoil - eair)/((Tsoil-ABSZERO)*Rtot);  // 0.0078 = 0.001*3600*18e-3/8.31 with 18e-3 = the molar mass of water vapor in kg/mol and 8.31 the ideal gas constant in J/mol/K; 0.001*3600*nbhours_covered*sites_per_dcell*LH*LH is used to convert evaporation in kg m-2 s-1 to m3 per day per dcell.
                
                //if (soil_phi3D[0][d] < -1) {
                // cout << "r_soil=" << r_soil << " r_soil_sellers=" << r_soil_sellers <<" r_aero=" <<r_aero << " VPDground=" << VPDground << " esat_ground=" << esat_ground << " esoil=" << esoil << " eair=" << eair << " soil_phi3D[0][d]=" << soil_phi3D[0][d] << " Tsoil=" << Tsoil << " TopWindSpeed_DCELL[d]=" << TopWindSpeed_DCELL[d] << " Wind ground level=" << exp(-alphaInoue*(1-1/Canopy_height_DCELL[d]))*TopWindSpeed_DCELL[d] << " evaporation S92=" << e  << " evaporation S92_sellers=" << e_sellers << " e_granier=" << e_Granier << " SWC3D[0][d]-Min_SWC[0]=" << SWC3D[0][d]-Min_SWC[0] << endl;
                // }
                
                Evaporation[d]=fmaxf(0.0,fminf(e,(SWC3D[0][d]-Min_SWC[0]))); // the amount of water evaporated from the soil cannot result in a water content below the residual water content. A model depending on soil matric potential would not need this.
                if(Evaporation[d]<0 || isnan(Evaporation[d])  || (SWC3D[0][d]-Min_SWC[0])<0 ) {
                    cout << "evaporation=" << Evaporation[d] << " and e=" << e << " and SWC3D[0][d]=" << SWC3D[0][d] << " and Min_SWC[0]=" << Min_SWC[0] << "and SWC3D[0][d]-Min_SWC[0]=" << SWC3D[0][d]-Min_SWC[0] << "; raero=" << r_aero << "; Canopy_height_DCELL[d]=" << Canopy_height_DCELL[d]<< "; TopWindSpeed_DCELL[d]=" << TopWindSpeed_DCELL[d] << endl;
                }
                
                SWC3D[0][d]-=Evaporation[d];
                
                // Refilling by rainfall
                
                Interception[d]=fminf(precip, 0.2*LAI_DCELL[0][d]);      // This is the amount of rainfall - in mm, as rainfall -, intercepted by vegetation cover, following the approach used in Liang et al. 1994 Journal of Geophysical Reserach, and also used by Laio et al. 2001 Advances in Water Resources and Fischer et al. 2014 Environmental Modelling & Software (FORMIX3, Madagascar). More complex approach can be used however - see eg. Gutierrez et al. 2014 Plos One (FORMIND, Chili), or Wagner et al. 2011 AFM (Paracou)
                Throughfall[d]=precip-Interception[d];
                Throughfall[d]*=sites_per_dcell*LH*LH*0.001; // to convert in absolute amount of water entering the soil voxel in m3
                
                if(isnan(Throughfall[d]) || (Throughfall[d]) <0 ) {
                    cout << "Incorrect throughfall" << endl;
                    cout << precip << "\t" << Interception[d] << "\t" << LAI_DCELL[0][d] << endl;
                }
                
                float in=Throughfall[d];
                
                /*if(SWC3D[0][d]<Max_SWC[0]) {
                    int l=0;
                    while((l<nblayers_soil) && (in>0.0)) {
                        if(in>(Max_SWC[l]-SWC3D[l][d])) {
                            in-=(Max_SWC[l]-SWC3D[l][d]);
                            SWC3D[l][d]=Max_SWC[l];
                            if(isnan(SWC3D[l][d]) || (SWC3D[l][d]-Min_SWC[l])<=0) {
                                cout << "incorrect SWC3D, Min/Max_SWC" << endl;
                                cout <<Max_SWC[l] << endl;
                            }
                        }
                        else{
                            SWC3D[l][d]+=in;
                            if (isnan(SWC3D[l][d]) || (SWC3D[l][d]-Min_SWC[l])<0) {
                                cout << "incorrect SWC3D, Min/Max_SWC" << endl;
                                cout << Throughfall[d] << "\t" <<in <<"\t" <<  precip << "\t" << Interception[d] << "\t" << LAI_DCELL[0][d] << endl;
                            }
                            in=0.0;
                        }
                        l++;
                    }
                }*/
                if(SWC3D[0][d]<Max_SWC[0]) {
                    int l=0;
                    while((l<nblayers_soil) && (in>0.0)) {
                        if(in>(FC_SWC[l]-SWC3D[l][d])) {
                            in-=(FC_SWC[l]-SWC3D[l][d]);
                            SWC3D[l][d]=FC_SWC[l];
                            if(isnan(SWC3D[l][d]) || (SWC3D[l][d]-Min_SWC[l])<=0) {
                                cout << "incorrect SWC3D, Min/Max_SWC" << endl;
                                cout <<Max_SWC[l] << endl;
                            }
                        }
                        else{
                            SWC3D[l][d]+=in;
                            if (isnan(SWC3D[l][d]) || (SWC3D[l][d]-Min_SWC[l])<0) {
                                cout << "incorrect SWC3D, Min/Max_SWC" << endl;
                                cout << Throughfall[d] << "\t" <<in <<"\t" <<  precip << "\t" << Interception[d] << "\t" << LAI_DCELL[0][d] << endl;
                            }
                            in=0.0;
                        }
                        l++;
                    }
                }
                else{ //if the top soil layer is already saturated (eg. inundated forest), throughfall -> runoff
                    Runoff[d]=Throughfall[d];
                }
                // Leakage
                Leakage[d]=in;
            }
            // END of the BUCKET MODEL.
            
            // Update of soil water potential field
            for (int d=0; d<nbdcells; d++) {
                for (int l=0; l<nblayers_soil; l++) {
                    //soil_phi3D[l][d]=phi_e[l]*pow((SWC3D[l][d]/Max_SWC[l]), -b[l]);
                    
                    float theta_w=(SWC3D[l][d]-Min_SWC[l])/(Max_SWC[l]-Min_SWC[l]);
                    
if (_WATER_RETENTION_CURVE==1) {
                    if(theta_w==0) {
                        theta_w=0.001; // SS addition for limit value
                        cout << "Warning theta_w = 0 " << endl ;
                    }
                    soil_phi3D[l][d]=a_vgm[l]*pow((pow(theta_w,-b_vgm[l])-1), c_vgm[l]); // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                    float inter= 1-pow((1-pow(theta_w, b_vgm[l])),m_vgm[l]);
                    Ks[l][d]=Ksat[l]*pow(theta_w, 0.5)*inter*inter; // this is the van Genuchten-Mualem model (as in Table 1 in Marthews et al. 2014)
                    
                    if (isnan(soil_phi3D[l][d]) || isnan(Ks[l][d]) ||  (SWC3D[l][d]-Min_SWC[l])<0) //|| KsPhi[l][d]==0.0 || Ks[l][d]==0.0 || soil_phi3D[l][d]==0.0)
                        cout << "In bucket model, layer " << l << " dcell " << d << " theta_w=" << theta_w << " SWC3D[l][d]-Min_SWC[l]=" << (SWC3D[l][d]-Min_SWC[l]) << " soil_phi3D[l][d]=" << soil_phi3D[l][d] << " Ksat=" << Ksat[l] << " Ks[l][d]=" << Ks[l][d] << endl ;
                    

} else if (_WATER_RETENTION_CURVE==0) {
                    soil_phi3D[l][d]=phi_e[l]*pow(theta_w, -b[l]); // this is the soil water characteristic of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                    Ks[l][d]=Ksat[l]*pow(theta_w, 2.5+2*b[l]); // this is the hydraulic conductivity curve of Brooks & Corey-Mualem (as in Table 1 in Marthews et al. 2014)
                    KsPhi[l][d]=Ksat[l]*phi_e[l]*pow(theta_w, 2.5+b[l]); //Ks times soil_phi3D, computed directly as the exact power of theta.
            
                    if (isnan(soil_phi3D[l][d]) || isnan(Ks[l][d]) ||  isnan(KsPhi[l][d]) || (SWC3D[l][d]-Min_SWC[l])<0) //|| KsPhi[l][d]==0.0 || Ks[l][d]==0.0 || soil_phi3D[l][d]==0.0)
                        cout << "In bucket model, layer " << l << " dcell " << d << " theta_w=" << theta_w << " SWC3D[l][d]-Min_SWC[l]=" << (SWC3D[l][d]-Min_SWC[l]) << " soil_phi3D[l][d]=" << soil_phi3D[l][d] << " Ksat=" << Ksat[l] << " phi_e=" << phi_e[l] <<" b[l]=" << b[l] << " KsPhi[l][d]=" << KsPhi[l][d] << " Ks[l][d]=" << Ks[l][d] << endl ;
                    //KsPhi2[l][d]=Ksat[l]*phi_e[l]*pow(theta_w, 2.5);
                    // we may want to shift to the van Genuchten-Mualem expressions of soil_phi3D and Ks, as the van genuchten-Mualem model is currently defacto the more standard soil hydraulic model (see ref in Table 1 in Marthews et al. 2014). To do so, see if we have data of soil pH, cation exchange capacity, organic carbon content, to explicitly compute the parameters with Hodnett & Tomasella 2002 (as recommended by Marthews et al. 2014 -- Table 2; or instead directly use the parameter provided by the map in Marthews et al. 2014.
                    
}
                }
            }
#endif
        }

            float layer_depth_previous = 0.0;
            for(int l=0; l<nblayers_soil; l++) {
                float soilWC=0.0;
                for (int d=0; d<nbdcells;d++) {
                    soilWC+=SWC3D[l][d]; // in m3
                }
                float layer_depth_current = layer_depth[l];
                float layer_thickness = layer_depth_current - layer_depth_previous;
                soilWC*=isites/layer_thickness;  // in m3/m3
                output[11] << soilWC << "\t";

                layer_depth_previous = layer_depth_current;
            }
            for(int l=0; l<nblayers_soil; l++) {
                float soilPhi=0.0;
                for(int d=0; d<nbdcells;d++) {
                    soilPhi+=soil_phi3D[l][d];  //in MPa
                }
                soilPhi*=icells; // in MPa
                output[11] << soilPhi << "\t";
            }
            
            output[11] <<"\n";

            if (iter==(nbiter-90) || iter==(nbiter-45) || iter==(nbiter-1)) {
                
                int o_swc, o_swp, o_wfluxes;
                if (iter==(nbiter-90)) {
                    o_swc=1;
                    o_swp=4;
                    o_wfluxes=18;
                }
                if (iter==(nbiter-45)) {
                    o_swc=2;
                    o_swp=5;
                    o_wfluxes=19;
                }
                if (iter==(nbiter-1)) {
                    o_swc=3;
                    o_swp=6;
                    o_wfluxes=20;
                }
                
                output[o_wfluxes] << "LAI" << "\t";
                for (int d=0; d<nbdcells;d++) {
                    output[o_wfluxes] << LAI_DCELL[0][d] << "\t";
                }
                output[o_wfluxes] << endl;
                
                output[o_wfluxes] << "Evaporation" << "\t";
                for (int d=0; d<nbdcells;d++) {
                    output[o_wfluxes] << Evaporation[d]*i_sites_per_dcell << "\t"; // in m
                }
                output[o_wfluxes] << endl;
                
                layer_depth_previous = 0.0;
                for(int l=0; l<nblayers_soil; l++) {
                    float layer_depth_current = layer_depth[l];
                    float layer_thickness = layer_depth_current - layer_depth_previous;
                    float norm=i_sites_per_dcell/layer_thickness;
                    output[o_swc] << l << "\t";
                    output[o_swp] << l << "\t";
                    output[o_wfluxes] << "Transpiration_" << l << "\t";
                    for (int d=0; d<nbdcells;d++) {
                        output[o_swc] << SWC3D[l][d]*norm  << "\t"; // in m3/m3
                        output[o_swp] << soil_phi3D[l][d] << "\t";
                        output[o_wfluxes] << Transpiration[l][d]*i_sites_per_dcell << "\t";
                    }
                    layer_depth_previous = layer_depth_current;
                    output[o_swc] << endl;
                    output[o_swp] << endl;
                    output[o_wfluxes] << endl;
                }
                
                
                
            }
            
            output[22] << iter  << "\t";
            output[23] << iter  << "\t";
            output[24] << iter  << "\t";
            for (int l=0; l<HEIGHT+1; l++) {
                LAI_young[l]*=isites;
                LAI_mature[l]*=isites;
                LAI_old[l]*=isites;
                output[22] << LAI_young[l]  << "\t";
                output[23] << LAI_mature[l]  << "\t";
                output[24] << LAI_old[l]  << "\t";
                LAI_young[l]=0.0;
                LAI_mature[l]=0.0;
                LAI_old[l]=0.0;
            }
            output[22] << endl;
            output[23] << endl;
            output[24] << endl;
            
            abund_phi_root*=inbhectares/sum1;
            abund10_phi_root*=inbhectares/sum10;
            agb_phi_root*=inbhectares/agb;
            
            output[31] << iter << "\t" << abund_phi_root << "\t" << abund10_phi_root << "\t" << agb_phi_root << endl;
            
            

            
#endif
