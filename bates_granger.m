file_path = '/Users/lakeshkansakar/aerosol_v2/all_final_qad_angles_cfrac.mat';
load(file_path);

number_sites = 50;
number_satellites = 5;
number_models = 31;

%locations = [2; 4; 6; 15; 20; 25; 50; 65; 80; 10; 11; 12; 13; 14; 15; 31; 35; 38; 71; 72];
locations = (1:100)';
years = [2005; 2006; 2007; 2008; 2009; 2010];

variance = zeros(number_sites,1);
covariance = zeros(10,1);

measurements_aeronet = [];

measurements_terra = [];
measurements_aqua = [];
measurements_misr = [];
measurements_omi = [];
measurements_seawifs = [];

total_days = [];

%Fetch the aeronet data and individual satellite data
for i=1:size(locations,1)
    for j=1:size(years,1)
        location_data = [];
        switch years(j)
            case 2005
                disp('2005');
                location_data = location_details_angles_qad_2005(locations(i,1));
                
            case 2006
                disp('2006');
                location_data = location_details_angles_qad_2005(locations(i,1));
                
            case 2007
                disp('2007');
                location_data = location_details_angles_qad_2005(locations(i,1));
                
            case 2008
                disp('2008');
                location_data = location_details_angles_qad_2005(locations(i,1));
                
            case 2009
                disp('2009');
                location_data = location_details_angles_qad_2005(locations(i,1));
                
            case 2010
                disp('2010');
                location_data = location_details_angles_qad_2005(locations(i,1));
        end
       
        total_days = [total_days; location_data.total_days_total];
        measurements_aeronet = [measurements_aeronet; location_data.measurements_aeronet];
        measurements_terra = [measurements_terra; location_data.measurements_terra];
        measurements_aqua = [measurements_aqua; location_data.measurements_aqua];
        measurements_misr = [measurements_misr; location_data.measurements_misr];
        measurements_omi = [measurements_omi; location_data.measurements_omi];
        measurements_seawifs = [measurements_seawifs; location_data.measurements_seawifs];
        
    end
end

measurements = [measurements_terra measurements_aqua measurements_misr measurements_omi measurements_seawifs];

%Calculate the error variance of individual satellites
for i=1:number_sites
    switch(i)
        case 1
            index = find(measurements_aeronet(:,:) > 0 & measurements_terra(:,:) > 0);
            available_measurements_satellite = measurements_terra(index,:);
            available_measurements_aeronet = measurements_aeronet(index,:);
            
        case 2
            index = find(measurements_aeronet(:,:) > 0 & measurements_aqua(:,:) > 0);
            available_measurements_satellite = measurements_aqua(index,:);
            available_measurements_aeronet = measurements_aeronet(index,:);
            
        case 3
            index = find(measurements_aeronet(:,:) > 0 & measurements_misr(:,:) > 0);
            available_measurements_satellite = measurements_misr(index,:);
            available_measurements_aeronet = measurements_aeronet(index,:);
            
        case 4
            index = find(measurements_aeronet(:,:) > 0 & measurements_omi(:,:) > 0);
            available_measurements_satellite = measurements_omi(index,:);
            available_measurements_aeronet = measurements_aeronet(index,:);
            
        case 5
            index = find(measurements_aeronet(:,:) > 0 & measurements_seawifs(:,:) > 0);
            available_measurements_satellite = measurements_seawifs(index,:);
            available_measurements_aeronet = measurements_aeronet(index,:);
    end
    
    %error = available_measurements_aeronet - available_measurements_satellite;
    %variance(i,1) = var(error);
    
end

%Calculate the covariance of errors between individual satellites
position = 1;
for i=1:number_satellites
    for j=1:number_satellites
        if(i ~= j && i< j)
            index = find(measurements(:,i) > 0 & measurements(:,j) > 0);
            data = measurements(index,:);
            covariance(position,1) = cov(data(:,i), data(:,j));
        end
    end
end

global satellite_data aeronet_data number_parameters;

for i=1:number_models
    switch(i)
        case 1
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:4) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            
        case 2
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 3 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            
        case 3
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:3) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 4
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1), (0.00000001));
            
        case 5
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 4]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 6
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 2 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 7
            index = find(measurements(:,1) <= 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1:2) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1 1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 8
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4 5]) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            
        case 9
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 4]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 10
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 11
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) <= 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 3]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1 1],(1),[ 0.00000001 0.00000001 0.00000001]);
            
        case 12
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 13
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1, 1],[],[],[1,1 1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 14
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[1 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 15
            index = find(measurements(:,1) <= 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,1) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            
        case 16
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:5) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            
        case 17
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:4) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 18
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 3 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 19
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) <= 0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2:3) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 20
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 21
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 22
            index = find(measurements(:,1) > 0 & measurements(:,2) <=0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[2 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);   
            
        case 23
            index = find(measurements(:,1) > 0 & measurements(:,2) <= 0 & measurements(:,3) >0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,2) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            
        case 24
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <=0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4 5]) = [];
            number_parameters = 2;
            k = fmincon(@find_parameters,[1,1],[],[],[1,1],(1),[0.00000001 0.00000001]);
            
        case 25
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) <= 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 4]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 26
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=0 & measurements(:,4) > 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,[3 5]) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 27
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) <=  0 & measurements(:,4) >0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,3) = [];
            number_parameters = 1;
            k = fmincon(@find_parameters,(1),[],[],(1),(1),(0.00000001));
            
        case 28
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <= 0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4:5) = [];
            number_parameters = 3;
            k = fmincon(@find_parameters,[1,1,1],[],[],[1,1,1],(1),[0.00000001 0.00000001 0.00000001]);
            
        case 29
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) <=0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,4) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            
        case 30
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) >0 & measurements(:,5) <= 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            satellite_data(:,5) = [];
            number_parameters = 4;
            k = fmincon(@find_parameters,[1,1,1,1],[],[],[1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001]);
            
        case 31
            index = find(measurements(:,1) > 0 & measurements(:,2) > 0 & measurements(:,3) > 0 & measurements(:,4) > 0 & measurements(:,5) > 0 & measurements_aeronet(:,1) > 0);
            satellite_data = measurements(index,:);
            aeronet_data = measurements_aeronet(index,:);
            number_parameters = 5;
            k = fmincon(@find_parameters,[1,1,1,1,1],[],[],[1,1,1,1,1],(1),[0.00000001 0.00000001 0.00000001 0.00000001 0.00000001]);
            
        
            
    end
end