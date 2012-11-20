%Implementation of the methods given in the Bates and Granger paper



%Load Nemanja's preprocessed data
file_path = '/Users/lakeshkansakar/aerosol_v2/all_final_qad_angles_cfrac.mat';
load(file_path);
load('location_details_aod_2005.mat');
load('location_details_aod_2006.mat');
load('location_details_aod_2007.mat');
load('location_details_aod_2008.mat');
%load('location_details_aod_2009.mat');
%load('location_details_aod_2010.mat');

%We will five satellite instruments(MISR, MODIS terra, MODIS aqua, OMI,
%Seawifs
number_satellites = 5;


%%%%%%%%%%%%%%%%% Method I (minimize the variance of the error)
%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
%%% Doesn't consider the weight assigned to the individual satellites to be
%%% the function of time i.e they are independent of time %%%%%%%%%%

%Use the 2^n-1 model where n is the number of satellites(instruments)
number_models = 31;

%Number of aeronet sites
number_locations = 171;


%years = [2005; 2006; 2007; 2008; 2009; 2010];
years = [2005; 2006; 2007; 2008];

%variance = zeros(number_sites,1);
%covariance = zeros(10,1);

measurements_aeronet = [];
measurements_terra = [];
measurements_aqua = [];
measurements_misr = [];
measurements_omi = [];
measurements_seawifs = [];

%Select just the sites at east and learn the model(test on the west side)
boundary_longitude = -900;

total_days = [];

%Fetch the aeronet data and individual satellite data
for i=1:number_locations
    for j=1:size(years,1)
        location_data = [];
        match_location_criteria = false;
        
        switch years(j)
            case 2005
                disp('2005');
                location_data = location_details_angles_qad_2005(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2005{1,i};
                    match_location_criteria = true;
                end
                
            case 2006
                disp('2006');
                location_data = location_details_angles_qad_2006(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2006{1,i};
                    match_location_criteria = true;
                end
                
            case 2007
                disp('2007');
                location_data = location_details_angles_qad_2007(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2007{1,i};
                    match_location_criteria = true;
                end
                
            case 2008
                disp('2008');
                location_data = location_details_angles_qad_2008(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2008{1,i};
                    match_location_criteria = true;
                end
                
            case 2009
                disp('2009');
                location_data = location_details_angles_qad_2009(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2009{1,i};
                    match_location_criteria = true;
                end
                
            case 2010
                disp('2010');
                location_data = location_details_angles_qad_2010(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    location_data = location_details_aod_2010{1,i};
                    match_location_criteria = true;
                end
        end
        
        if(match_location_criteria == true && ~isempty(location_data.measurements))
            try 
                
                measurements_aeronet = [measurements_aeronet; location_data.measurements(1,:)'];
                measurements_terra = [measurements_terra; location_data.measurements(2,:)'];
                measurements_misr = [measurements_misr; location_data.measurements(3,:)'];
                measurements_aqua = [measurements_aqua; location_data.measurements(4,:)'];
                measurements_omi = [measurements_omi; location_data.measurements(5,:)'];
                measurements_seawifs = [measurements_seawifs; location_data.measurements(6,:)'];

                match_location_criteria = false;
            catch exception
                pause();
            end
        end
        
    end
    
end

measurements = [measurements_terra measurements_aqua measurements_misr measurements_omi measurements_seawifs];

global satellite_data aeronet_data number_parameters;

%Train the models based upon the satellite data combination(2^n-1 combinations)
%TODO - find an elegant way to do this
%Use one year out cross validation for testing

parameters = zeros(2^number_satellites-1,number_satellites+1);

for i=1:number_models
    switch(i)
        case 1
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:4) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            parameters(i,number_satellites) = k;
            parameters(i,6) = size(index,1);
            
        case 2
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 3 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            parameters(i,number_satellites-1) = k;
            parameters(i,6) = size(index,1);
            
        case 3
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:3) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,4:5) = k;
            parameters(i,6) = size(index,1);
            
        case 4
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1), (0.00000001));
            parameters(i,3) = k;
            parameters(i,6) = size(index,1);
            
        case 5
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,3) = k(1);
            parameters(i,5) = k(2);
            parameters(i,6) = size(index,1);
            
        case 6
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,3) = k(1);
            parameters(i,4) = k(2);
            parameters(i,6) = size(index,1);
            
        case 7
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:2) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1 1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,3) = k(1);
            parameters(i,4) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 8
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            parameters(i,2) = k;
            parameters(i,6) = size(index,1);
            
        case 9
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,5) = k(2);
            parameters(i,6) = size(index,1);
            
        case 10
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,4) = k(2);
            parameters(i,6) = size(index,1);
            
        case 11
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1 1],(1),[ 0.00000001 0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,4) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 12
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,3) = k(2);
            parameters(i,6) = size(index,1);
            
        case 13
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1, 1],[],[],[1,1 1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,3) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 14
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,3) = k(2);
            parameters(i,4) = k(3);
            parameters(i,6) = size(index,1);
            
        case 15
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,2) = k(1);
            parameters(i,3) = k(2);
            parameters(i,4) = k(3);
            parameters(i,5) = k(4);
            parameters(i,6) = size(index,1);
            
        case 16
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:5) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            parameters(i,1) = k;
            parameters(i,6) = size(index,1);
            
        case 17
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:4) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,5) = k(2);
            parameters(i,6) = size(index,1);
            
        case 18
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 3 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,4) = k(2);
            parameters(i,6) = size(index,1);
            
        case 19
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <= 0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:3) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,4) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 20
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,3) = k(2);
            parameters(i,6) = size(index,1);
            
        case 21
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,3) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 22
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);   
            parameters(i,1) = k(1);
            parameters(i,3) = k(2);
            parameters(i,4) = k(3);
            parameters(i,6) = size(index,1);
            
        case 23
            index = find(measurements(:,1) > 0 & measurements(:,2) <= 0 & measurements(:,3) >0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,3) = k(2);
            parameters(i,4) = k(3);
            parameters(i,5) = k(4);
            parameters(i,6) = size(index,1);
            
            
        case 24
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,6) = size(index,1);
            
        case 25
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,5) = k(3);
            parameters(i,6) = size(index,1);
            
        case 26
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,4) = k(3);
            parameters(i,6) = size(index,1);
            
        case 27
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=  0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,3) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,4) = k(3);
            parameters(i,5) = k(4);
            parameters(i,6) = size(index,1);
            
        case 28
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4:5) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,3) = k(3);
            parameters(i,6) = size(index,1);
            
        case 29
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,3) = k(3);
            parameters(i,5) = k(4);
            parameters(i,6) = size(index,1);
            
        case 30
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) >0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,5) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,3) = k(3);
            parameters(i,4) = k(4);
            parameters(i,6) = size(index,1);
            
        case 31
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            number_parameters = 5;
            k = fmincon(@find_parameters,[1,1,1,1,1],[],[],[1,1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001 0.00000001]);
            parameters(i,1) = k(1);
            parameters(i,2) = k(2);
            parameters(i,3) = k(3);
            parameters(i,4) = k(4);
            parameters(i,5) = k(5);
            parameters(i,6) = size(index,1);
    end
end




















%%%%%%%%%%%%%%%%%%%%%%% Method II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Consider the varying weight feature i.e the weight depends upon
%%%%%%%%% time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%TODO - There is lots of redundancy in this code. Clean this up by one
%%function

%Select just the sites at east and learn the model(test on the west side)
boundary_longitude = -88;

total_days = [];

%Window size to calculate the error from. The weight at timestamp T depends upon the mistakes 
%made by the predictor in the previous 10 timestamps
window_size = 10;
years = [2005; 2006; 2007; 2008; 2009; 2010];

loc_data = {};
measurements_aeronet_all_sites = [];
measurements_terra_all_sites = [];
measurements_aqua_all_sites = [];
measurements_misr_all_sites = [];
measurements_omi_all_sites = [];
measurements_seawifs_all_sites = [];
location_aod_all_years_instruments = {};

for i=1:number_locations
    for j=1:size(years,1)
        location_data = [];
        match_location_criteria = false;

        switch years(j)
            case 2005
                disp('2005');
                location_data = location_details_angles_qad_2005(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2005.mat');
                end

            case 2006
                disp('2006');
                location_data = location_details_angles_qad_2006(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2006.mat');
                end

            case 2007
                disp('2007');
                location_data = location_details_angles_qad_2007(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2007.mat');
                end

            case 2008
                disp('2008');
                location_data = location_details_angles_qad_2008(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2008.mat');
                end

            case 2009
                disp('2009');
                location_data = location_details_angles_qad_2009(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2009.mat');
                end

            case 2010
                disp('2010');
                location_data = location_details_angles_qad_2010(1,i);
                location = location_data.location;
                if(location(1,1) > boundary_longitude)
                    match_location_criteria = true;
                    location_compact_aod = load('location_details_aod_2010.mat');
                end
        end

        if(match_location_criteria == true && ~isempty(location_compact_aod.measurements))
            total_days_total = location_data.total_days_total;
            start_day_total = location_data.start_day_total;
            end_day_total = location_data.end_day_total;
            day_of_year = dayofyear(start_day_total(1),start_day_total(2),start_day_total(3));
            measurements_aeronet = zeros(366,1);
            measurements_terra = zeros(366,1);
            measurements_aqua = zeros(366,1);
            measurements_misr = zeros(366,1);
            measurements_omi = zeros(366,1);
            measurements_seawifs = zeros(366,1);
            measurements_aeronet(day_of_year:total_days_total) = location_compact_aod.measurements(:,1)';
            measurements_terra(day_of_year:total_days_total) = location_compact_aod.measurements(:,2)';
            measurements_aqua(day_of_year:total_days_total) = location_compact_aod.measurements(:,3)';
            measurements_misr(day_of_year:total_days_total) = location_compact_aod.measurements(:,4)';
            measurements_omi(day_of_year:total_days_total) = location_compact_aod.measurements(:,5)';
            measurements_seawifs(day_of_year:total_days_total) = location_compact_aod.measurements(:,6)';
            
            match_location_criteria = false;
        end
        
        measurements_aeronet_all_years = [measurements_aeronet_all_years; measurements_aeronet];
        measurements_terra_all_years = [measurements_terra_all_years; measurements_terra];
        measurements_aqua_all_years = [measurements_aqua_all_years; measurements_aqua];
        measurements_misr_all_years = [measurements_misr_all_years; measurements_misr];
        measurements_omi_all_years = [measurements_omi_all_years; measurements_omi];
        measurements_seawifs_all_years = [measurements_seawifs_all_years; measurements_seawifs];
        
        location_aod_all_years_instruments{i}.measurements_aeronet_all_years = measurements_aeronet_all_years;
        location_aod_all_years_instruments{i}.measurements_terra_all_years = measurements_terra_all_years;
        location_aod_all_years_instruments{i}.measurements_aqua_all_years = measurements_aqua_all_years;
        location_aod_all_years_instruments{i}.measurements_misr_all_years = measurements_misr_all_years;
        location_aod_all_years_instruments{i}.measurements_omi_all_years = measurements_omi_all_years;
        location_aod_all_years_instruments{i}.measurements_seawifs_all_years = measurements_seawifs_all_years;
    end

end


%%%%%%%%%%%%%%%%%%%%%%% Method III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
