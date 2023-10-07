%% %%%%%%%%%%%%%%%%%%%% %%
%% %%%%% HDGM_sim %%%%% %%
%% %%%%%%%%%%%%%%%%%%%% %%

%%% Simulates a HDGM model

function [DSTEM_obj_sim,DSTEM_obj_sim_stem_par,...
    DSTEM_str_ground_sim, DSTEM_ground] = DSTEM_HDGM_sim(DSTEM_str_ground,...
    DSTEM_str_sim_setup,n_covs,n_sites,X_names,X_varcov,...
    datestamp_begin, datestamp_end, spatialCV_step)

d1 = datetime(datestamp_begin,'InputFormat','dd-MM-yyyy HH:mm');
d2 = datetime(datestamp_end,'InputFormat','dd-MM-yyyy HH:mm');
time_stamps = d1:days(1):d2;
T = length(time_stamps);
    
%%% Computation begin
date_begin = datetime('now')


%% %%%%% HDGM data building
poll = DSTEM_str_ground.poll;
n = cell(length(poll),1);
sites = cell(length(poll),1);
obj_stem_gridlist_p = stem_gridlist();
for p = 1:length(poll)
    %%% Sampling sites
    idx_st_poll = ismember(DSTEM_str_ground.ARPA_stats_reg.New_cod_stz,...
        DSTEM_str_ground.(['Coords_' DSTEM_str_ground.poll{p}]).IDStat);
    SOI = find(idx_st_poll);
    sites{p} = datasample(SOI,n_sites{1,p},'Replace',false);
    sites_idx{p} = find(ismember(SOI,sites{p}));    
    %%% Dependent variables
    ground.Y{p} = DSTEM_str_ground.([poll{p}])(sites_idx{p},1:T);
    ground.Y_name{p} = poll{p};
    n{p} = n_sites{1,p};
    %%% Loading covariates for the selected pollutants at each monitoring stations
    ground.X_beta_name{p} = ['Constant' X_names];
    ground.X_beta{p} = zeros(n{p},n_covs+1,T);
    for tt = 1:T
        ground.X_beta{p}(:,:,tt) = [repmat(1,n{p},1) , ...
            mvnrnd(zeros(n_covs,1),X_varcov{p},n{p})];
		% ground.X_beta{p}(:,:,tt) = [repmat(1,n{p},1) , ...
        %     mvnrnd(ones(n_covs,1),X_varcov{p},n{p})];
    end
    %%% Coordinates grid
    ground.coordinates{p} = [DSTEM_str_ground.(['Coords_' poll{p}]).Latitude(sites_idx{p}), ...
        DSTEM_str_ground.(['Coords_' poll{p}]).Longitude(sites_idx{p})];
    %%% Constant term
    ground.X_z{p} = ones(n{p}, 1);
    ground.X_z_name{p} = {['constant_' poll{p}]};
    %%% STEM_grid
    % contains all the information related to the sampling locations of
    % a single variable
    obj_stem_grid = cell(length(poll),1);
    obj_stem_grid{p} = stem_grid(ground.coordinates{p}, 'deg', 'sparse', 'point');
    %%% STEM_gridlist
    % collector of the stem_grid objects for all the variables;
    obj_stem_gridlist_p.add(obj_stem_grid{p});
    if spatialCV_step ~= 0
        %%% Cross-validation settings
        S_val{1,p} = 1:spatialCV_step:n{p};
        ground.S_val{1,p} = S_val{1,p};
    end
end


%% %%%%%
ground.obj_stem_gridlist_p = obj_stem_gridlist_p;
ground.datestamp_begin = datestamp_begin;
ground.datestamp_end = datestamp_end;
DSTEM_ground = ground;


%% %%%%% STEM settings
%%% STEM_varset
% Contains the observed data of all the variables and the loading coefficients;
obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);
%%% datestamp
% contains the information related to the date and time of the observations;
% obj_stem_datestamp = stem_datestamp('01-01-2017 00:00','31-12-2017 00:00',T);
obj_stem_datestamp = stem_datestamp(datestamp_begin,datestamp_end,T);
%%% Shapefile
shape = [];
%%% Model type
obj_stem_modeltype = stem_modeltype('HDGM');


stem_data_input.stem_modeltype = obj_stem_modeltype;
stem_data_input.stem_varset_p = obj_stem_varset_p;
stem_data_input.stem_gridlist_p = obj_stem_gridlist_p;
stem_data_input.stem_datestamp = obj_stem_datestamp;
stem_data_input.shape = shape;
if spatialCV_step ~= 0
   obj_stem_validation = stem_validation(ground.Y_name,S_val,0,...
       repmat({'point'},1,length(poll)));
   stem_data_input.stem_validation = obj_stem_validation;
end
obj_stem_data = stem_data(stem_data_input);


% if spatialCV_step ~= 0
%    obj_stem_validation = stem_validation(ground.Y_name,S_val,0,...
%        repmat({'point'},1,length(poll)));
%    obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
%                              [], [], obj_stem_datestamp, obj_stem_validation, ...
%                              obj_stem_modeltype, shape);
% else
%    obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
%        [], [], obj_stem_datestamp, [], ...
%        obj_stem_modeltype, shape);
% end


%%% DSTEM_par object creation
% Contains the structure and the values of the model parameters;
obj_stem_par_constraints = stem_par_constraints();
%%% time_diagonal
% Used to specify if the matrices G and Sigma_eta are diagonal or not
obj_stem_par_constraints.time_diagonal = 0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
%%% DSTEM_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
%%% DSTEM_par export
DSTEM_obj_sim_stem_par = obj_stem_par;


%% %%%%% Simulation
DSTEM_obj_sim = stem_sim(obj_stem_model);
% Questo va messo prima, altrimenti incorrelato sempre...
DSTEM_obj_sim.stem_model.stem_par.beta = DSTEM_str_sim_setup.beta;
DSTEM_obj_sim.nan_rate = DSTEM_str_sim_setup.nan_rate;
DSTEM_obj_sim.nan_pattern_par = DSTEM_str_sim_setup.nan_pattern_par;

DSTEM_obj_sim.stem_model.stem_par_initial = DSTEM_obj_sim.stem_model.stem_par;
DSTEM_obj_sim.stem_model.stem_par_initial.theta_z = DSTEM_str_sim_setup.theta_z;
DSTEM_obj_sim.stem_model.stem_par_initial.v_z = DSTEM_str_sim_setup.v_z;
DSTEM_obj_sim.stem_model.stem_par_initial.sigma_eta = DSTEM_str_sim_setup.sigma_eta;
DSTEM_obj_sim.stem_model.stem_par_initial.G = DSTEM_str_sim_setup.G;
DSTEM_obj_sim.stem_model.stem_par_initial.sigma_eps = DSTEM_str_sim_setup.sigma_eps;
simulate(DSTEM_obj_sim,DSTEM_obj_sim.nan_rate,DSTEM_obj_sim.nan_pattern_par)


DSTEM_str_ground_sim = DSTEM_str_ground;
% contains(fieldnames(DSTEM_str_ground_sim),'X_')
fnames = fieldnames(DSTEM_str_ground_sim);
DSTEM_str_ground_sim = rmfield(DSTEM_str_ground_sim,fnames(contains(fieldnames(DSTEM_str_ground_sim),{'NO2','NOx','PM10','PM2_5'})));
% fnames = fieldnames(DSTEM_str_ground_sim);
% DSTEM_str_ground_sim = rmfield(DSTEM_str_ground_sim,fnames(contains(fieldnames(DSTEM_str_ground_sim),'Coords_')));
DSTEM_str_ground_sim.ARPA_stats_reg = DSTEM_str_ground.ARPA_stats_reg(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.lon = DSTEM_str_ground.lon(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.lat = DSTEM_str_ground.lat(unique(vertcat(sites{:,1})),:);
DSTEM_str_ground_sim.sites = size(DSTEM_str_ground_sim.ARPA_stats_reg,1);
DSTEM_str_ground_sim.var_names = DSTEM_str_ground.poll;
DSTEM_str_ground_sim.measure_type = 'simulated';
d1 = datenum(datetime(datestamp_begin,'InputFormat','dd-MM-yyyy HH:mm'));
d2 = datenum(datetime(datestamp_end,'InputFormat','dd-MM-yyyy HH:mm'));
DSTEM_str_ground_sim.date_time = d1:d2;
DSTEM_str_ground_sim.time_stamps = T;
for p = 1:length(poll)
    DSTEM_str_ground_sim.(['Coords_' poll{p}]) = DSTEM_str_ground.(['Coords_' poll{p}])(sites_idx{1,p},:);
    DSTEM_str_ground_sim.(poll{p}) = DSTEM_obj_sim.stem_model.stem_data.stem_varset_p.Y{p,1};
    DSTEM_str_ground_sim.(['X_' poll{p}]) = DSTEM_obj_sim.stem_model.stem_data.stem_varset_p.X_beta{1,p};
end


return;

end