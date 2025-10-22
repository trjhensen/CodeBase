function OK = runMarsADRC(outputPathMars, marsInputTablePath)
% Run mars using tracktable configuration

% Prerequisites: 
% This step is done by performing the following script
% Go to root directory in the terminal or powershell
% # Download repository
% git clone https://www.github.com/ThieleLAB/mars-pipeline
% 
% # Go to directory
% cd mars-pipeline
% 
% # In a debian based linux: first do:
% python3 -m venv .venv
% source .venv/bin/activate
% python3 -m pip install -r requirements.txt
% 
% # Install requirements
% pip install -r requirements.txt
%
% OPTIONAL:
% # Run WSL and do the rest in wsl
% # wsl --install # install ubuntu
% # wsl
% 
% # build streamlit app in docker
% sudo docker build -t my_streamlit_app .
% 
% # Run docker app
% sudo docker run -p 8501:8501 my_streamlit_app
% 
% # Open app in browser
% http://localhost:8501
% # Upload MARS input files and run in MARS


OK = false; % Preallocate output

% Force delete all current outputs in the processedMicrobiomePath
%if isfolder(outputPathMars); rmdir(outputPathMars,'s'); end

% Character array variable to the folder where the output of MARS is stored. 
if ~isfolder(outputPathMars)
    mkdir(outputPathMars)
end


% *REQUIRED*. Character array variable with path to the microbiome taxonomy 
% and read abundance file. 

% *REQUIRED*. Character array variable to the directory where the 
% mars-pipeline is present on the system.
paths.Mars.marsRepoPath = 'C:\Users\mspg\Documents\mars-pipeline';

% *REQUIRED*. Character array variable to the file that starts python.
paths.Mars.pythonPath = 'C:\Users\mspg\anaconda3\python.exe';

% *REQUIRED*. Character array variable with path to the microbiome taxonomy 
% and read abundance file.
paths.Mars.readsTablePath = marsInputTablePath;

% Character array variable indicating the desired file format for saving 
% outputs: this is needed for the relative abundance file path used in the
% next line.
paths.Mars.outputExtensionMars = 'csv';

% Path to relative abundance file. Defaults to MARS output, update if not using MARS.
paths.Mars.relAbunFilePath = fullfile(outputPathMars,'renormalized_mapped_forModelling',['renormalized_mapped_forModelling_species.', paths.Mars.outputExtensionMars]);

% Numeric value for total read counts per sample under which samples are
% excluded from analysis. Only applies when readsTable contains absolute
% read counts (not relative abundances). Defaults to 1, with a minimum of 1.
paths.Mars.sample_read_counts_cutoff = 1;

% Numeric value under which relative abundances are considered to be zero.
paths.Mars.cutoffMars = 0.000001;

% String to the file where OTUs are matched to taxonomic assignments.
% OPTIONAL if the taxonomic assignments are already in the readsTable. 
% REQUIRED if not.
paths.Mars.taxaTable = string(missing);

% A boolean to indicate if the genus name is in the name of the species e.g.
% Prevotella copri. If genus name is in species name, set to false. 
% Otherwise, set to true. OPTIONAL, defaults to false.
paths.Mars.flagLoneSpecies = false;

% The delimiter used to separate taxonomic levels.
paths.Mars.taxaDelimiter = ';';

% A boolean specifying if one wants to remove clade name extensions from
% all taxonomic levels of microbiome taxa. If set to false, MARS might find
% significantly fewer models in AGORA2, as clade extensions are not included there.
paths.Mars.removeClade = true;

% A string defining if AGORA2, APOLLO, a combination of both, or a user-defined
% database should be used as the model database to check presence in. 
% Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db", "user_db".
paths.Mars.reconstructionDb = "full_db";

% A string containing the full path to the user-defined database,
% which should be in .csv, .txt, .parquet, or .xlsx format and
% have column names = taxonomic levels. Only required if whichModelDatabase
% is set to "user_db".
paths.Mars.userDbPath = "";

disp(' > Perform metagenomic mapping using MARS.');
disp(' > Initialise python environment.');

% Perform metagenomic mapping
runMars(paths.Mars.pythonPath, ...
    paths.Mars.marsRepoPath, ...
    paths.Mars.cutoffMars, ...
    paths.Mars.outputExtensionMars, ...
    paths.Mars.flagLoneSpecies, ...
    paths.Mars.taxaDelimiter, ...
    paths.Mars.removeClade, ...
    paths.Mars.reconstructionDb, ...
    paths.Mars.userDbPath, ...
    paths.Mars.sample_read_counts_cutoff, ...
    paths.Mars.readsTablePath, ...
    paths.Mars.taxaTable, ...
    outputPathMars)

OK = true;
end