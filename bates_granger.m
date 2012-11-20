%Implementation of the methods given in the Bates and Granger paper



%Load Nemanja's preprocessed data
file_path = '/Users/lakeshkansakar/aerosol_v2/all_final_qad_angles_cfrac.mat';
load(file_path);
load('location_details_aod_2005.mat');
load('location_details_aod_2006.mat');
load('location_details_aod_2007.mat');
load('location_details_aod_2008.mat');
load('location_details_aod_2009.mat');
load('location_details_aod_2010.mat');

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


years = [2005; 2006; 2007; 2008; 2009; 2010];

%variance = zeros(number_sites,1);
%covariance = zeros(10,1);

measurements_aeronet = [];
measurements_terra = [];
measurements_aqua = [];
measurements_misr = [];
measurements_omi = [];
measurements_seawifs = [];

%Select just the sites at east and learn the model(test on the west side)
boundary_longitude = -88;

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

% %Calculate the error variance of individual satellites
% for i=1:number_sites
%     switch(i)
%         case 1
%             index = find(measurements_aeronet(:,:) > 0 & measurements_terra(:,:) > 0);
%             available_measurements_satellite = measurements_terra(index,:);
%             available_measurements_aeronet = measurements_aeronet(index,:);
%             
%         case 2
%             index = find(measurements_aeronet(:,:) > 0 & measurements_aqua(:,:) > 0);
%             available_measurements_satellite = measurements_aqua(index,:);
%             available_measurements_aeronet = measurements_aeronet(index,:);
%             
%         case 3
%             index = find(measurements_aeronet(:,:) > 0 & measurements_misr(:,:) > 0);
%             available_measurements_satellite = measurements_misr(index,:);
%             available_measurements_aeronet = measurements_aeronet(index,:);
%             
%         case 4
%             index = find(measurements_aeronet(:,:) > 0 & measurements_omi(:,:) > 0);
%             available_measurements_satellite = measurements_omi(index,:);
%             available_measurements_aeronet = measurements_aeronet(index,:);
%             
%         case 5
%             index = find(measurements_aeronet(:,:) > 0 & measurements_seawifs(:,:) > 0);
%             available_measurements_satellite = measurements_seawifs(index,:);
%             available_measurements_aeronet = measurements_aeronet(index,:);
%     end
%     
%     %error = available_measurements_aeronet - available_measurements_satellite;
%     %variance(i,1) = var(error);
%     
% end

%Calculate the covariance of errors between individual satellites
% position = 1;
% for i=1:number_satellites
%     for j=1:number_satellites
%         if(i ~= j && i< j)
%             index = find(measurements(:,i) > 0 & measurements(:,j) > 0);
%             data = measurements(index,:);
%             covariance(position,1) = cov(data(:,i), data(:,j));
%         end
%     end
% end

global satellite_data aeronet_data number_parameters;

%Train the models based upon the satellite data combination(2^n-1 combinations)
%TODO - find an elegant way to do this
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
%%%%%%%%% time


measurements_aeronet = [];
measurements_terra = [];
measurements_aqua = [];
measurements_misr = [];
measurements_omi = [];
measurements_seawifs = [];

%Select just the sites at east and learn the model(test on the west side)
boundary_longitude = -88;

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
    
    measurements_aeronet_all_sites{i} = measurements_aeronet;
    measurements_terra_all_sites{i} = measurements_terra;
    measurements_misr_all_sites{i} = measurements_misr;
    measurements_aqua_all_sites{i} = measurements_aqua;
    measurements_omi_all_sites{i} = measurements_omi;
    measurements_seawifs_all_sites{i} = measurements_seawifs;
    
end


%%TO DO - There is lots of redundancy in this code. Clean this up by one
%%function
parameters = zeros(2^number_satellites-1,number_satellites+1);
window = 5;

for i=1:number_models
    switch(i)
        case 1
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:4) = [];
            number_parameters = 1;
            parameters{i}.parameters = ones(size(aeronet_data(:,1),1));
            parameters{i}.number_points = size(index,1);
            
        case 2
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 3 5]) = [];
            number_parameters = 1;
            parameters{i}.parameters = ones(size(aeronet_data(:,1),1));
            parameters{i}.number_points = size(index,1);
            
        case 3
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:3) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 4
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4 5]) = [];
            number_parameters = 1;
            parameters{i}.parameters = ones(size(aeronet_data(:,1),1));
            parameters{i}.number_points = size(index,1);
            
        case 5
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 6
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 7
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:2) = [];
            number_parameters = 3;
            
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 8
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4 5]) = [];
            number_parameters = 1;
            parameters{i}.parameters = ones(size(aeronet_data(:,1),1));
            parameters{i}.number_points = size(index,1);
            
        case 9
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 10
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 11
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 12
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 13
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 14
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 5]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 15
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1) = [];
            number_parameters = 4;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            
            k1(1,1) = 1/4;
            k2(1,1) = 1/4;
            k3(1,1) = 1/4;
            k4(1,1) = 1/4;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4)/(E1+E2+E3+E4);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4)/(E1+E2+E3+E4);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4)/(E1+E2+E3+E4);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3)/(E1+E2+E3+E4);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4];
            parameters{i}.number_points = size(index,1);
            
        case 16
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:5) = [];
            number_parameters = 1;
            parameters{i}.parameters = ones(size(aeronet_data(:,1),1));
            parameters{i}.number_points = size(index,1);
            
        case 17
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:4) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 18
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 3 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 19
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <= 0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:3) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 20
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 21
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 22
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 5]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 23
            index = find(measurements(:,1) > 0 & measurements(:,2) <= 0 & measurements(:,3) >0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2) = [];
            number_parameters = 4;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            
            k1(1,1) = 1/4;
            k2(1,1) = 1/4;
            k3(1,1) = 1/4;
            k4(1,1) = 1/4;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4)/(E1+E2+E3+E4);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4)/(E1+E2+E3+E4);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4)/(E1+E2+E3+E4);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3)/(E1+E2+E3+E4);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4];
            parameters{i}.number_points = size(index,1);
            
            
        case 24
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4 5]) = [];
            number_parameters = 2;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k1(1,1) = 0.5;
            k2(1,1) = k1(1,1);
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E2 = sum(error2_windowed(j:j-window));
                E1 = sum(error1_windowed(j:j-window));
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * E2/(E1+E2);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * E1/(E1+E2);
            end
                
            parameters{i}.parameters = [k1 k2];
            parameters{i}.number_points = size(index,1);
            
        case 25
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 26
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 5]) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 27
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=  0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,3) = [];
            number_parameters = 4;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            
            k1(1,1) = 1/4;
            k2(1,1) = 1/4;
            k3(1,1) = 1/4;
            k4(1,1) = 1/4;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4)/(E1+E2+E3+E4);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4)/(E1+E2+E3+E4);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4)/(E1+E2+E3+E4);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3)/(E1+E2+E3+E4);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4];
            parameters{i}.number_points = size(index,1);
            
        case 28
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4:5) = [];
            number_parameters = 3;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            
            k1(1,1) = 1/3;
            k2(1,1) = 1/3;
            k3(1,1) = 1/3;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3)/(E1+E2+E3);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3)/(E1+E2+E3);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2)/(E1+E2+E3);
            end
                
            parameters{i}.parameters = [k1 k2 k3];
            parameters{i}.number_points = size(index,1);
            
        case 29
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4) = [];
            number_parameters = 4;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            
            k1(1,1) = 1/4;
            k2(1,1) = 1/4;
            k3(1,1) = 1/4;
            k4(1,1) = 1/4;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4)/(E1+E2+E3+E4);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4)/(E1+E2+E3+E4);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4)/(E1+E2+E3+E4);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3)/(E1+E2+E3+E4);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4];
            parameters{i}.number_points = size(index,1);
            
        case 30
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) >0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,5) = [];
            number_parameters = 4;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            
            k1(1,1) = 1/4;
            k2(1,1) = 1/4;
            k3(1,1) = 1/4;
            k4(1,1) = 1/4;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4)/(E1+E2+E3+E4);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4)/(E1+E2+E3+E4);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4)/(E1+E2+E3+E4);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3)/(E1+E2+E3+E4);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4];
            parameters{i}.number_points = size(index,1);
            
        case 31
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            number_parameters = 5;
            error1 = abs(aeronet_data - satellite_data(:,1));
            error2 = abs(aeronet_data - satellite_data(:,2));
            error3 = abs(aeronet_data - satellite_data(:,3));
            error4 = abs(aeronet_data - satellite_data(:,4));
            error5 = abs(aeronet_data - satellite_data(:,5));
            
            a = 1;
            b = ones(1,4);
            error1_windowed = filter(b,a,error1);
            error2_windowed = filter(b,a,error2);
            error3_windowed = filter(b,a,error3);
            error4_windowed = filter(b,a,error4);
            error5_windowed = filter(b,a,error5);
            
            k1 = zeros(size(error1,1),1);
            k2 = k1;
            k3 = k1;
            k4 = k1;
            k5 = k1;
            
            k1(1,1) = 1/5;
            k2(1,1) = 1/5;
            k3(1,1) = 1/5;
            k4(1,1) = 1/5;
            k5(1,1) = 1/5;
            
            lambda = 0.5;
            
            for j=2:size(error1,1)
                E1 = sum(error1_windowed(j:j-window));
                E2 = sum(error2_windowed(j:j-window));
                E3 = sum(error3_windowed(j:j-window));
                E4 = sum(error4_windowed(j:j-window));
                E5 = sum(error5_windowed(j:j-window));
                
                k1(j,1) = lambda*k1(j-1) + (1-lambda) * (E2+E3+E4+E5)/(E1+E2+E3+E4+E5);
                k2(j,1) = lambda*k2(j-1) + (1-lambda) * (E1+E3+E4+E5)/(E1+E2+E3+E4+E5);
                k3(j,1) = lambda*k3(j-1) + (1-lambda) * (E1+E2+E4+E5)/(E1+E2+E3+E4+E5);
                k4(j,1) = lambda*k4(j-1) + (1-lambda) * (E1+E2+E3+E5)/(E1+E2+E3+E4+E5);
                k4(j,1) = lambda*k5(j-1) + (1-lambda) * (E1+E2+E3+E4)/(E1+E2+E3+E4+E5);
            end
                
            parameters{i}.parameters = [k1 k2 k3 k4 k5];
            parameters{i}.number_points = size(index,1);
    end
end



%%%%%%%%%%%%%%%%%%%%%%% Method III %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
