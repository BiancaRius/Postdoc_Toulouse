# 🌿 TROLL Model Documentation


## Main function

The main function serves as the entry point of the TROLL model. **It initializes** the model, processes input parameters, and executes the main simulation loop.

### Overview and Description 

#### **1. Initial Setup and Configuration**  

##### **MPI Initialization (if enabled)**  
- Checks if the code was compiled with MPI support (`#ifdef MPI`).  
- If enabled:  
  - Initializes the MPI environment (`MPI_Init`).  
  - Retrieves the process rank (`MPI_Comm_rank`) and the total number of processes (`MPI_Comm_size`).  
- If MPI is not enabled:  
  - Sets `mpi_rank = 0` and `mpi_size = 1`, indicating single-process execution.  

##### **Command-Line Argument Processing**  
- Iterates through the program’s command-line arguments (`argc`, `argv`).  
- Identifies and extracts input and output file names based on specific flags:  
  - `-i`: General input file.  
  - `-o`: Output file.  
  - `-m`: Climate input file.  
  - `-s`: Species input file.  
  - `-f`: Inventory data file.  
- Stores file names in corresponding variables (`inputfile`, `inputfile_climate`, `outputfile`, etc.).  

##### **Input File Reading and Parameter Initialization**  
- Reads configuration data and model parameters from the input files using:  
  - `ReadInputGeneral()`  
  - `ReadInputPointcloud()`  
  - `ReadInputInventory()`  
- Initializes the random number generator (`gsl_rng`).  
- Creates a log file for simulation details.  
- Calls initialization functions:  
  - `Initialise()` – Sets up global model parameters and call input functions.
  - `InitialiseOutputStreams()` – Prepares output streams.  
  - `AllocMem()` – Allocates memory for dynamic variables.  
- If the ABC module is enabled, it initializes related data.  

##### **Initial Information Display**  
- Prints simulation parameters such as CO2 concentration, species count, and active modules.  

---

#### **2. Simulation Evolution Loop**  

###### **Main Simulation Loop**  
- Runs a `for` loop for `nbiter` iterations, representing simulation time steps.  
- In each iteration:  
  - Calls `Evolution()` to update the model state.  
  - Measures execution time for performance analysis.  
  - If visualization is enabled, calls `OutputVisual()` at predefined intervals.  
  - If the ABC module is enabled, updates transmittance and outputs using:  
    - `UpdateTransmittanceCHM_ABC()`  
    - `OutputABC()`  
    - `UpdateDBHtrackingABC()`  

##### **Intermediate and Final Data Output**  
- Generates simulation state snapshots (`OutputSnapshot()`) at key time points (initial, intermediate, final).  
- If extended outputs are enabled, calls:  
  - `OutputLAI()` – Outputs Leaf Area Index data.  
  - `OutputCHM()` – Outputs Canopy Height Model data.  
- If inventory output is enabled, saves soil water content data.  

---

#### **3. Finalization and Cleanup**  

##### **Final Information Display**  
- Prints the number of surviving trees at the end of the simulation.  
- Calculates and displays total execution time.  

##### **Memory Cleanup and Program Termination**  
- Closes output files (`CloseOutputs()`).  
- Releases allocated memory (`FreeMem()`).  
- If MPI is enabled, finalizes the MPI environment (`MPI_Finalize()`).  
- Terminates the program (`exit()`).  

---

### **Function Call Breakdown**  

#### **Initialization and Configuration Functions**  
- `MPI_Init()`: Initializes the MPI environment (if enabled).  
- `MPI_Comm_rank()`: Retrieves the MPI process rank.  
- `MPI_Comm_size()`: Retrieves the number of MPI processes.  
- `atoi()`: Converts strings to integers (used in argument parsing).  
- `sprintf()`: Formats and stores strings (used for filenames).  
- `ReadInputGeneral()`: Reads general parameters from the input file.  
- `gsl_rng_env_setup()`: Sets up the GSL random number generator environment.  
- `gsl_rng_default()`: Retrieves the default GSL RNG type.  
- `gsl_rng_alloc()`: Allocates memory for the GSL RNG.  
- `gsl_rng_set()`: Sets the RNG seed.  
- `Initialise()`: Initializes global model parameters.  
- `InitialiseOutputStreams()`: Prepares output streams.  
- `AllocMem()`: Allocates memory for model variables.  
- `InitialiseABC()`: Initializes the ABC module (if enabled).  
- `ReadInputPointcloud()`: Reads parameters for point cloud generation.  
- `ReadInputInventory()`: Reads initial inventory data.  

#### **Evolution Loop Functions**  
- `clock()`: Measures execution time.  
- `Evolution()`: Advances the model by one time step.  
- `GetTimeofyear()`: Retrieves the current simulation time within the year.  
- `OutputVisual()`: Generates visual outputs.  
- `ExportPointcloud()`: Exports point cloud data (if enabled).  
- `UpdateTransmittanceCHM_ABC()`: Updates transmittance for the ABC module.  
- `OutputABC()`: Outputs ABC module data.  
- `UpdateDBHtrackingABC()`: Updates DBH tracking in the ABC module.  

#### **Data Output Functions**  
- `OutputSnapshot()`: Creates simulation snapshots at key points.  
- `OutputLAI()`: Outputs Leaf Area Index data.  
- `OutputCHM()`: Outputs Canopy Height Model data.  

#### **Finalization Functions**  
- `fmaxf()`: Returns the maximum of two floating-point values.  
- `MPI_Reduce()`: Performs MPI reduction operations (if enabled).  
- `CloseOutputs()`: Closes output files.  
- `FreeMem()`: Frees allocated memory.  
- `MPI_Finalize()`: Finalizes the MPI environment.  
- `exit()`: Terminates the program.  


### OBS:
```cpp
for(int argn=1; argn<argc; argn++){ // Arguments of the input and output files
    if(*argv[argn] == '-'){
        switch(*(argv[argn]+1)){
Explanation of Each Component:
for(int argn=1; argn<argc; argn++):
```

This is a for loop that iterates through the command-line arguments passed to the program.
argn is initialized to 1 because argv[0] is typically the name of the program itself, and we want to skip it to process the actual arguments.
The loop continues as long as argn is less than argc, which is the total number of command-line arguments.
After each iteration, argn is incremented by 1.
if(*argv[argn] == '-'):

This condition checks if the current argument (pointed to by argv[argn]) starts with a '-'.
In command-line arguments, options or flags usually start with a '-', indicating that this argument is a command-line option rather than a regular input or output file.
switch(*(argv[argn]+1)):

This is a switch statement that evaluates the character immediately following the '-' in the command-line argument.
*(argv[argn] + 1) dereferences the pointer to get the second character of the argument string (the one right after the '-').
The switch statement allows you to handle different command-line options based on what follows the '-'.


### Inputs

The input files, besides general parameter are called inside function `Initialise()`, for example: `ReadInputSpecies()`, `ReadInputClimate()`, `ReadInputDailyvar` and `ReadInputSoil`

#### Soil inputs

Function used `ReadInputSoil()`

##### Description
The ReadInputSoil function reads soil property data from an input file and processes it for further use in the model. It extracts information on soil texture (sand, silt, and clay proportions), organic carbon content, bulk density, pH, and cation exchange capacity (CEC), among other parameters. The function also computes key soil properties, such as saturated soil water content (Sat_SWC) and residual soil water content (Res_SWC), using pedotransfer functions based on the selected water retention curve model.

The function processes soil properties layer by layer and applies empirical formulas to estimate soil hydraulic properties.
If additional field data (e.g., pH, organic content, bulk density, CEC) is available, it can be used to refine soil parameter estimations.
If no local soil data is provided, default soil parameters should be taken from global datasets, as suggested in Marthews et al. (2014).
The function logs the number of soil layers read and prints computed values for verification.

##### Functionalities
- Opens and reads the soil input file (inputfile_soil).

- Extracts the following soil layer attributes:
    - Layer thickness (m) 
    - Texture proportions: Silt (%), Clay (%), Sand (%)
    - Soil organic carbon (SOC) (% → converted to gC/kg)
    - Dry bulk density (DBD) (g/cm³)
    - pH
    - Cation exchange capacity (CEC) (cmol/kg)

- Computes cumulative soil depth for each layer.

- Calculates soil hydrological properties based on the _WATER_RETENTION_CURVE setting:
    If _WATER_RETENTION_CURVE == 1: Uses Hodnett & Tomasella (2002) functions for tropical soils.
    If _WATER_RETENTION_CURVE == 0: Uses Tomasella & Hodnett (1998) texture-based pedotransfer functions.

- Computes:
    - Saturated soil water content (Sat_SWC) (m³/m³)
    - Maximum soil water capacity (Max_SWC) (m³)
    - Field capacity (FC_SWC)
    - Saturated hydraulic conductivity (Ksat) (mm/s or kg/m²/s)
    - Residual soil water content (Res_SWC)
    - Minimum soil water capacity (Min_SWC) (m³)

- Allocates memory dynamically for soil parameter arrays and checks for allocation errors.

##### Function (code and description)

```cpp
        void ReadInputSoil(){
            cout << endl << "Reading in file: " << inputfile_soil << endl;
            
            /// @brief Open the soil input file for reading
            fstream InSoil(inputfile_soil, ios::in);
            /// fstream is a library to manage files both for reading and writing
            /// it creates an objet fstream type that is called InSoil
            /// inputfile_soil is the name of the file that is being open
            /// ios::in specifies it will be open in the input mode (for reading)
            
            if(InSoil){
                InSoil.getline(buffer,256,'\n');
                /// InSoil.getline is a member of function fstream class used for reading a line of text
                /// buffer: is a character array and 256 is the maximum of characters that will be read
                /// '\n': getline() will read characters from the file until it encounters a newline character 

                // Soil property vectors
                /// declare vector that will store the values for the variables
                vector<float> layer_thickness, proportion_Silt, proportion_Clay, proportion_Sand; // in m, %, %,%
                vector<float> SOC, DBD, pH, CEC; //soil organic content, provided in %; dry bulk density, in g cm-3; pH; cation exchange capacity, in cmol kg-1
                
                // Reserve space for vectors
                /// at leats 20 spaces for storing
                SOC.reserve(20);
                DBD.reserve(20);
                pH.reserve(20);
                CEC.reserve(20);
                layer_thickness.reserve(20);
                proportion_Silt.reserve(20);
                proportion_Clay.reserve(20);
                proportion_Sand.reserve(20);

                // Initializes the variable that will count the number of soil layers
                nblayers_soil = 0;
                
                // we go through all lines in the input file
                string line;
                /// getline(InSoil, line):reads a line of text from the InSoil file stream and stores it in the line variable.
                while(getline(InSoil, line)){
                    /// Process each line:
                    istringstream linestream(line);

// Read different sets of variables depending on _WATER_RETENTION_CURVE flag                    
if (_WATER_RETENTION_CURVE==1) {
                    /// declare variables of float type (they store temporary values that is read in a line of the file)
                    float thickness_current, proportion_Silt_current, proportion_Clay_current, proportion_Sand_current, SOC_current, DBD_current, pH_current, CEC_current;

                    /// linestream is an istringstream object that cotains the line of data read from file
                    /// >> extract values from linestream and converts it to float before storing it
                    /// the reading occurs in the order the variables are listed
                    linestream >> thickness_current >> proportion_Silt_current >> proportion_Clay_current >> proportion_Sand_current >> SOC_current >> DBD_current >> pH_current >> CEC_current;
                    
                    /// @brief Store the read values into respective vectors
                    /// push_back() adds a new element at the end of the vector
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

                    /// @brief Store the read values into respective vectors
                    layer_thickness.push_back(thickness_current);
                    proportion_Silt.push_back(proportion_Silt_current);
                    proportion_Clay.push_back(proportion_Clay_current);
                    proportion_Sand.push_back(proportion_Sand_current);
}
                    /// when the loop ends (it starts at 236) it means that all file lines were read, then
                    /// nblayers_soil will containd the number of read lines that corresponds to the total number
                    /// of soil layers found in the file
                    nblayers_soil++;
                }

                cout << "Read in: " << nblayers_soil << " soil layers" << endl;
                
                // Memory allocation for the var layer_depth and depth calculation
                if(NULL==(layer_depth=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc layer_depth" << endl;
                float cumulative_depth=0.0;
                
                /// @brief Compute cumulative depth for each soil layer
                for (int l=0; l<nblayers_soil; l++) {
                    cumulative_depth+=layer_thickness[l];
                    layer_depth[l]=cumulative_depth;
                }

                // Compute soil water characteristics
                // (added in the header but kept in here) in this version, all soil parameters (Sat_SWC, Res_SWC) are computed from soil texture data (%clay, %silt, %sand) provided in input for each layer. If additional information is available from the field (soil pH, organic content, dry bulk density, cation exchange capacity), this could be also provided in input and used to refine the computation of these soil parameters (see Table 2 in Marthews et al. 2014 Geoscientific Model Development and Hodnett & Tomasella 2002 Geoderma -- for tropical soils). Alternatively, if no local field soil data is available, these soil parameters (Sat_SWC, Res_SWC) should be drawn from global maps and databases --see Marthews et al. 2014, and directly provided in input. ==> ccl: to standardize the input file, the soil parameters (Sat_SWC, Res_SWC) should probably be provided in input, and the computation of those properties from the available local data (here %clay, %silt, %sand) made using a new function of RconTROLL (and not here)
                // (added in the header but kept in here) Sat_SWC and Res_SWC are here computed according Tomasella & Hodnett 1998 from soil texture information (see Table 2 in Marthews et al. 2014)
                if(NULL==(Sat_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Sat_SW" << endl;
                if(NULL==(Max_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Max_SW" << endl;
                if(NULL==(FC_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc FC_SW" << endl;

                /// for each soil layer:
                for (int l=0; l<nblayers_soil; l++) {
                    
if (_WATER_RETENTION_CURVE==1) {
                    /// @brief Compute saturation soil water content using Hodnett & Tomasella 2002 formula
                    /// Water content in saturation (water volume/soil volume) - the fraction os soil volume that is occupied by water
                    Sat_SWC[l] = 0.01*(81.799+(0.099*proportion_Clay[l])-(31.42*DBD[l])+(0.018*CEC[l])+(0.451*pH[l])-(0.0005*proportion_Sand[l]*proportion_Clay[l]));  // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3

} else if (_WATER_RETENTION_CURVE==0) {
                    /// @brief Compute saturation soil water content using Tomasella & Hodnett 1998 formula
                    Sat_SWC[l]= 0.01*(40.61+(0.165*proportion_Silt[l])+(0.162*proportion_Clay[l])+(0.00137*proportion_Silt[l]*proportion_Silt[l])+(0.000018*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3
}
                    /// The maximum water content (total amount of water)
                    Max_SWC[l]=Sat_SWC[l]*sites_per_dcell*LH*LH*layer_thickness[l]; // in m3

                    cout << "layer " << l << " Vol=" << sites_per_dcell*LH*LH*layer_thickness[l]<< " m3; Sat_SWC =" << Sat_SWC[l] << " MAX_SWC =" << Max_SWC[l] << " m3." << endl;
                }
                
                // Compute saturated hydraulic conductivity
                if(NULL==(Ksat=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Ksat" << endl;
                for (int l=0; l<nblayers_soil; l++) {
                    /// @brief Compute hydraulic conductivity using Cosby et al. 1984 equation
                    Ksat[l]=0.007055556*pow(10,(-0.60-(0.0064*proportion_Clay[l])+(0.0126*proportion_Sand[l]))); // according to Cosby et al. 1984 (the only expression of k_sat reported in Table 2 of Marthews et al. 2014). k_sat is here in mm/s or equivalently in kg/m2/s.
                    cout << "layer " << l << " Ksat=" << Ksat[l] <<"mm/s or kg/m2/s  "<< Ksat[l]*9.8/18 << endl;
                }
                
                // Compute residual soil water content
                if(NULL==(Res_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Res_SW" << endl;
                if(NULL==(Min_SWC=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc Max_SW" << endl;
                
                // Loop through each soil layer
                for (int l=0; l<nblayers_soil; l++) {
                    
if (_WATER_RETENTION_CURVE==1) {
                    /// @brief Compute residual soil water content using Hodnett & Tomasella 2002 formula
                    /// minimum amount of water that stays in the soil even when the majority of available water
                    /// was already removed by evaporation, plant transpiration and runoff                    
                    Res_SWC[l]= 0.01*(22.733-(0.164*proportion_Sand[l])+(0.235*CEC[l])-(0.831*pH[l])+(0.0018*proportion_Clay[l]*proportion_Clay[l])+(0.0026*proportion_Sand[l]*proportion_Clay[l])); // this is the Hodnett & Tomasella 2002 tropical pedotransfer function, as reported in Table 2 of Marthews et al. 2014. in m3.m-3

} if (_WATER_RETENTION_CURVE==0) {
                    /// @brief Compute residual soil water content using Tomasella & Hodnett 1998 formula
                    Res_SWC[l]= 0.01*fmaxf(0.0,(-2.094+(0.047*proportion_Silt[l])+(0.431*proportion_Clay[l])-(0.00827*proportion_Silt[l]*proportion_Clay[l]))); // this is the Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
}
                    /// @brief Computes minimum soil water content in cubic meters
                    Min_SWC[l]=Res_SWC[l]*sites_per_dcell*LH*LH*layer_thickness[l]; //in m3

                    cout << "layer " << l << " Vol=" << sites_per_dcell*LH*LH*layer_thickness[l] << "m3; Res=" << Res_SWC[l]<<  " MIN_SWC =" << Min_SWC[l] << " m3" << endl;
                }
                
if (_WATER_RETENTION_CURVE==1) {
                if(NULL==(a_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc a_vgm" << endl;
                if(NULL==(b_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc b_vgm" << endl;
                if(NULL==(c_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc c_vgm" << endl;
                if(NULL==(m_vgm=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc m_vgm" << endl;

                
                for (int l=0; l<nblayers_soil; l++) {
                    /** 
                     * @brief Computes the alpha parameter for the van Genuchten-Mualem model 
                     * @details This calculation follows the Hodnett & Tomasella (2002) approach,
                     * as referenced in Table 2 of Marthews et al. (2014).
                     */
                    float alpha=1000*exp((-2.294-(3.526*proportion_Silt[l])+(2.440*(0.1*SOC[l]))-(0.076*CEC[l])-(11.331*pH[l])+(0.019*proportion_Silt[l]*proportion_Silt[l]))*0.01); // this is the alpha parameter of the van Genuchten-Mualem model, in MPa-1 (ie after already dividing by rho*g, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                    a_vgm[l]=-1.0/alpha;
                    
                    /** @brief Computes the n, m, b, and c parameters for the van Genuchten-Mualem model */
                    float n_vgm=exp((62.986-(0.833*proportion_Clay[l])-(0.529*(SOC[l]*0.1))+(0.593*pH[l])+(0.007*proportion_Clay[l]*proportion_Clay[l])-(0.014*proportion_Sand[l]*proportion_Silt[l]))*0.01);    // this is the n parameter of the van Genuchten-Mualem model, dimensionless, following Hodnett & Tomasella 2002, see Table 2 in Marthews et al. 2014
                    m_vgm[l]=1.0-1.0/n_vgm;
                    b_vgm[l]=1.0/m_vgm[l];
                    c_vgm[l]=1.0-m_vgm[l];
                    
                    /** 
                     * @brief Computes soil water content at field capacity (FC) - the amount of water retained in the soil
                     * after the excess of water is drained due to gravity. In the field capacity, soil water retention (matricial force)
                     * equilibrates the gravity force
                     * The water that stays in the soil is available to plants, once it is not as strongly retained as the residual water
                     * @details Uses a standard approach from Marthews et al. (2014),
                     * considering FC at -10 kPa instead of -33 kPa.
                     */
                    FC_SWC[l]=(Res_SWC[l] + (Sat_SWC[l]-Res_SWC[l])*pow((pow(0.01*alpha, 1/c_vgm[l])+1),-(1/b_vgm[l])))*sites_per_dcell*LH*LH*layer_thickness[l]; // this is the layer water content at field capacity, in m3. As in Marthews et al. 2014 (cf. note in Table 2), we used Phi at FC=-10kPa and not -33kPa, following Marshall et al., 1996; Townend et al., 2001; Tomasella and Hodnett, 2004)
                    
                    cout << "layer " << l << " alpha=" << alpha << "\t" << "n_vgm=" << n_vgm << " FC_SWC=" << FC_SWC[l] << endl;
                }

/** @brief Computes water retention parameters using Tomasella & Hodnett (1998) formula */
} else if (_WATER_RETENTION_CURVE==0) {
                if(NULL==(phi_e=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc phi_e" << endl;
                if(NULL==(b=new float[nblayers_soil])) cerr<<"!!! Mem_Alloc b" << endl;
                
                for (int l=0; l<nblayers_soil; l++) {
                    /**
                     * @brief Computes air-entry potential (phi_e) using Tomasella & Hodnett (1998)
                     * @details Formula derived from Marthews et al. (2014), based on tropical soils.
                     */
                    phi_e[l]=-0.001*(0.285+(0.000733*proportion_Silt[l]*proportion_Silt[l])-(0.00013*proportion_Silt[l]*proportion_Clay[l])+(0.0000036*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014. phi_e is here in MPa.
                    
                    /** @brief Computes parameter b for the water retention curve */
                    b[l]=exp(1.197+(0.00417*proportion_Silt[l])-(0.0045*proportion_Clay[l])+(0.000894*proportion_Silt[l]*proportion_Clay[l])-(0.00001*proportion_Silt[l]*proportion_Silt[l]*proportion_Clay[l])); // according to Tomasella & Hodnett 1998 tropical texture-based pedotransfer function, as reported in Table 2 of Marthews et al. 2014.
                    
                    //phi_e[l]=-0.00000001*pow(10.0,(2.17-(0.0063*proportion_Clay[l])-(0.0158*proportion_Sand[l])))*(1000*9.80665); // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. In MPa.
                    //b[l]=3.10+0.157*proportion_Clay[l]-0.003*proportion_Sand[l]; // according to Cosby et al. 1984, non -tropical and texture-based but widely used, as reported in Table 2 in Marthews et al. 2014. Dimensionless.
                    
                    /** 
                     * @brief Computes soil water content at field capacity (FC)
                     * @details Uses Marthews et al. (2014) approach, considering FC at -10 kPa.
                     */
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
```
### Water availability

Function used `Water_availability()`

#### Description

### Update Field (bucket model, update seeds, recruit Tree)
Function used `UpdateField()`

#### Description
The UpdateField() function advances the simulation by one time step. It updates daily climate variables, models neighborhood-based ecological interactions, handles tree recruitment, calculates the field's leaf area index (LAI), and simulates soil water dynamics (if enabled), including processes like interception, throughfall, evaporation, transpiration, runoff, and leakage. It also updates soil water potential based on the soil water content and a chosen water retention model.

#### Functionalities
- Its primary responsibility is to update various environmental and biological state variables within the simulated field for a single time step.

- Climate Data Assimilation:
    - reads daily climate variables (night temperature, precipitation, daily mean wind speed, daily mean irradiance, daily mean temperature, and daily mean vapor pressure deficit) from arrays.
    - The modulo operator (%) is used to cycle through the daily data.

- Possibility of declaring reduced precipitation (and others) experiments

- call UpdateSeeds()

- Neighborhood Density Dependence (NDD):

      If _NDD is true, it calculates a neighborhood density dependence field (t_NDDfield) for each site and species. It iterates through each site and then through a defined neighborhood (radius Rndd) around that site. This code block calculates, for each location in the field and for each species, a measure of neighborhood influence. This influence is based on the basal area of neighbors within a specific radius. The idea is that the presence and size of neighboring plants can affect the growth, survival, or other processes occurring at the central location.
      For each "site" in the field, it iterates through its neighbors within a defined radius (`Rndd`). If a neighbor is established (has a positive age), its basal area (proportional to the square of its diameter, `t_dbh`) contributes to the `t_NDDfield` of the central site, specific to the neighbor's species (`t_sp_lab`). This contribution is scaled by a normalization factor (`normBA`). Essentially, it quantifies how the size and proximity of neighboring entities influence each location in the field for each species present.


- Recruitment:

      It calls a function RecruitTree(), which handles the introduction of new individuals (e.g., seedlings) into the simulated field.

- Leaf Area Index (LAI) Calculation:

      It initializes a 3D leaf area index field (LAI3D).
      It calls a CalcLAI() method for each tree (T[site].CalcLAI()), implying that individual trees contribute to the overall LAI at different heights.
      It then cumulatively sums the LAI from the canopy top downwards, resulting in a vertical profile of LAI.

- Water Cycle Processes (if WATER is defined):

      It initializes various hydrological variables for each "dcell": runoff, interception, throughfall, evaporation, leakage, and  it also initializes LAI/canopy height related variables for dcells.

- LAI(average) at height h in dcell
      It calculates the average LAI at different heights (LAI_DCELL) and the average canopy   height (Canopy_height_DCELL) for each dcell based on the trees within that dcell.
      It estimates the wind speed at the top of the canopy (TopWindSpeed_DCELL) based on the daily mean wind speed and the canopy height, using a logarithmic or exponential relationship depending on whether the canopy height is below or above a meteorological station height.

- Soil water dynamics using a "bucket model":
      - Water Uptake: Reduces soil water content (SWC3D) based on tree transpiration (Transpiration).

      - Evaporation: Calculates and subtracts evaporation from the top soil layer based on a physically based model involving vapor pressure deficits, soil and air temperatures, and resistances (soil and aerodynamic).

      - Refilling by Rainfall: Calculates interception by the canopy and adds the remaining throughfall to the soil water content, prioritizing upper soil layers. If the topsoil is saturated, excess water becomes runoff.

      - Leakage: Accounts for water leaving the bottom of the soil profile.

      - Water Table Effect: If _WATER_TABLE is enabled, it sets the soil water content in layers below the water table depth (WTD) to the maximum water holding capacity.

      - It updates the soil water potential (soil_phi3D) and hydraulic conductivity (Ks, KsPhi) based on the current soil water content


#### Function (code and description)
```cpp
        //#################################
        // Global function: Update all fields
        //#################################
        //! - This is an important function for TROLL -- Includes many of the operations
        //! - set the iteration environment -- nb: the current structure of code suppose that environment is periodic (a period = a year), if one wants to input a variable climate, with interannual variation and climate change along the simulation, a full climatic input needs to be input (ie number of columns=iter and not iterperyear) and change iterperyear by nbiter here.
        void UpdateField() {
          #ifdef FULL_CLIMATE

    /**
    * @brief Assigns daily climate variables.
    *
    * This code snippet assigns values from daily climate data arrays to individual variables.
    * The modulo operator (%) is used to cycle through the daily data.
    *
    * @param iter      The current iteration number.
    * @param nbdays    The total number of days in the climate data cycle.
    */            
            tnight=NightTemperature[iter%nbdays];
            precip=Rainfall[iter%nbdays];
            WSDailyMean=DailyMeanWindSpeed[iter%nbdays];
            WDailyMean=DailyMeanIrradiance[iter%nbdays]*SWtoPPFD;
            tDailyMean=DailyMeanTemperature[iter%nbdays];
            VPDDailyMean=DailyMeanVapourPressureDeficit[iter%nbdays];

            // cout << "regular prec :  " << precip << endl;
            // Applying reduced precipitation for tests with WTD
            precip *= 0.3;
            //cout << "reduced prec :  " << precip << endl;
           
#else
            tnight=NightTemperature[iter%iterperyear];
            precip=Rainfall[iter%iterperyear];
            WSDailyMean=DailyMeanWindSpeed[iter%iterperyear];
            WDailyMean=DailyMeanIrradiance[iter%iterperyear]*SWtoPPFD;
            tDailyMean=DailyMeanTemperature[iter%iterperyear];
            VPDDailyMean=DailyMeanVapourPressureDeficit[iter%iterperyear];
            
#endif // FULL_CLIMATE
....... COMPLETE
        }
```
#### **Functions called**

##### 1. UpdateSeeds()
######  *1.a. Description*
######  *1.b. Functionalities*
######  *1.c. Function (code and description)*
######  *1.d. Functions called*

###### **1.d.1. DisperseSeed**
######  1.d.1.a. Description
######  1.d.1.b. Functionalities
######  1.d.1.c. Function (code and description)


##### 2. RecreuitTree()
######  *2.a. Description*
######  *2.b. Functionalities*
######  *2.c. Function (code and description)*

##### 3. CalcLAI()
######  *3.a. Description*
This is not a function but it is an important part of the code, so I decided to put it here as a function description
######  *3.b. Functionalities*
######  *3.c. Function (code and description)*

##### 4. Bucket model (NOT a function)
######  *4.a. Description*
This is not a function but it is an important part of the code, so I decided to put it here as a function description
######  *4.b. Functionalities*
######  *4.c. Function (code and description)*
#### Bucket model

for dcells

......... COMPLETE

calculates water uptake from transpiration
discount from SWC3D
calculates evaporation from soil
discount from SWCD3

refill the soil with rainfall (the amount of rainfall that goes to the soil)
calculates interception based on LAI
throughfall is converted to m3
soil layer weight

rename the variable
```cpp
in = Throughfall[d] 
```

##### **\_Water table depth - fixed values version\_**

This version introduces the option to include a fixed water table depth (WTD) in the model, with three configuration modes:

1. No WTD – original soil water dynamics.

2. Shallow WTD – last three soil layers are saturated.

3. Deep WTD – only the last soil layer is saturated.

- Possible tests:
  - vary layers depths
  - vary layers soil composition

###### Water Table Depth Definitions
Soil layers from surface to bottom:

| Layer | Thickness (m) | Cumulative Depth (m) |
|-------|----------------|----------------------|
| 1     | 0.10           | 0.10                 |
| 2     | 0.23           | 0.33                 |
| 3     | 0.40           | 0.73                 |
| 4     | 0.80           | 1.53                 |
| 5     | 0.97           | 2.50                 |

Water table depth settings:

- Shallow WTD: below the top two layers → WTD = 0.33 m

- Deep WTD: below the top four layers → WTD = 1.53 m

![Texto alternativo](comparingWTD_impementation/wtd_representation.png)


###### Code implementations

1. WTD on and off 

Control the water table feature via the a parameter in the input_global file:

```cpp
_WATER_TABLE = 0 // disables water table
_WATER_TABLE = 1 // enables water table
```
If _WATER_TABLE = 1, the model will simulate soil water content based on a fixed water table depth.

2. Soil water saturation depending on WTD

- Inside bucket model:

```cpp
if (_WATER_TABLE == 1){ /// WTD on
  int l=0; // layer counter
  while((l<nblayers_soil)) {
    //if the depth of the layer is higher than the wtd, the amount of water in the layer (SWC3D) is = max of water the soil layer can hold
    if(layer_depth[l]>WTD){
      SWC3D[l][d] = Max_SWC[l];
    }
      l++;
  }
}

```

**3. Comparing soil water content with and without WTD implementation**

![Texto alternativo](comparingWTD_impementation/hill.png)








 


