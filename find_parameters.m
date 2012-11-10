function [f df] = find_parameters(params)
  
    
    global satellite_data aeronet_data number_parameters;

    switch(number_parameters)
        case 1
            error = satellite_data - aeronet_data;
            k1 = params(1);
            f = k1^2 * var(error);
            df = 2*k1*var(error);
            
            
        case 2
            error1 = satellite_data(:,1) - aeronet_data;
            error2 = satellite_data(:,2) - aeronet_data;
            
            k1 = params(1);
            k2 = params(2);
            
            covariance = cov(error1, error2);
            f = k1^2*var(error1) + k2^2*var(error2) + 2*k1*k2*covariance(1,2);
            
            derivative_k1 = 2*k1*var(error1) + 2*k2*covariance(1,2);
            derivative_k2 = 2*k2*var(error2) + 2*k1*covariance(1,2);
            df = [derivative_k1 derivative_k2]; 
            
        case 3
            error1 = satellite_data(:,1) - aeronet_data;
            error2 = satellite_data(:,2) - aeronet_data;
            error3 = satellite_data(:,3) - aeronet_data;
            
            k1 = params(1);
            k2 = params(2);
            k3 = params(3);
            
            covariance12 = cov(error1, error2);
            covariance13 = cov(error1, error3);
            covariance23 = cov(error2, error3);
            
            f = k1^2*var(error1) + k2^2*var(error2) + k3^2*var(error3) + 2*k1*k2*covariance12(1,2) + 2*k1*k3*covariance13(1,3) + 2*k2*k3*covariance23(1,2);
            
            derivative_k1 = 2*k1*var(error1) + 2*k2*covariance12(1,2) + 2*k3*covariance13(1,2) ;
            derivative_k2 = 2*k2*var(error2) + 2*k1*covariance12(1,2) + 2*k3*covariance23(1,2) ;
            derivative_k3 = 2*k3*var(error3) + 2*k2*covariance23(1,2) + 2*k1*covariance13(1,2) ;
            df = [derivative_k1 derivative_k2 derivative_k3];
            
        case 4
            error1 = satellite_data(:,1) - aeronet_data;
            error2 = satellite_data(:,2) - aeronet_data;
            error3 = satellite_data(:,3) - aeronet_data;
            error4 = satellite_data(:,4) - aeronet_data;
            
            k1 = params(1);
            k2 = params(2);
            k3 = params(3);
            k4 = params(4);
            
            covariance12 = cov(error1, error2);
            covariance13 = cov(error1, error3);
            covariance14 = cov(error1, error4);
            covariance23 = cov(error2, error3);
            covariance24 = cov(error2, error4);
            covariance34 = cov(error3, error4);
            
            f = k1^2*var(error1) + k2^2*var(error2) + k3^2*var(error3) + k4^2*var(error4) + 2*k1*k2*covariance12(1,2) + 2*k1*k3*covariance13(1,3) + 2*k1*k4*covariance14(1,2) + ...
                2*k2*k3*covariance23(1,2) + 2*k2*k4*covariance24(1,2) + 2*k3*k4*covariance34(1,2);
            
            derivative_k1 = 2*k1*var(error1) + 2*k2*covariance12(1,2) + 2*k3*covariance13(1,2) + 2*k4*covariance14(1,2);
            derivative_k2 = 2*k2*var(error2) + 2*k1*covariance12(1,2) + 2*k3*covariance23(1,2) + 2*k4*covariance24(1,2);
            derivative_k3 = 2*k3*var(error3) + 2*k2*covariance23(1,2) + 2*k1*covariance13(1,2) + 2*k4*covariance34(1,2);
            derivative_k4 = 2*k4*var(error4) + 2*k2*covariance24(1,2) + 2*k1*covariance14(1,2) + 2*k3*covariance34(1,2);
            
            df = [derivative_k1 derivative_k2 derivative_k3 derivative_k4];
            
        case 5
            error1 = satellite_data(:,1) - aeronet_data;
            error2 = satellite_data(:,2) - aeronet_data;
            error3 = satellite_data(:,3) - aeronet_data;
            error4 = satellite_data(:,4) - aeronet_data;
            error5 = satellite_data(:,5) - aeronet_data;
            
            k1 = params(1);
            k2 = params(2);
            k3 = params(3);
            k4 = params(4);
            k5 = params(5);
            
            covariance12 = cov(error1, error2);
            covariance13 = cov(error1, error3);
            covariance14 = cov(error1, error4);
            covariance15 = cov(error1, error5);
            covariance23 = cov(error2, error3);
            covariance24 = cov(error2, error4);
            covariance25 = cov(error2, error5);
            covariance34 = cov(error3, error4);
            covariance35 = cov(error3, error5);
            covariance45 = cov(error4, error5);
            
            f = k1^2*var(error1) + k2^2*var(error2) + k3^2*var(error3) + k4^2*var(error4) + k5^2*var(error5) + 2*k1*k2*covariance12(1,2) + 2*k1*k3*covariance13(1,3) + 2*k1*k4*covariance14(1,2) + ...
                + 2*k1*k5*covariance15(1,2) +  2*k2*k3*covariance23(1,2) + 2*k2*k4*covariance24(1,2) + 2*k2*k5*covariance25(1,2) + 2*k3*k4*covariance34(1,2) + 2*k3*k5*covariance35(1,2) + ...
                + 2*k4*k5*covariance45(1,2);
            
            derivative_k1 = 2*k1*var(error1) + 2*k2*covariance12(1,2) + 2*k3*covariance13(1,2) + 2*k4*covariance14(1,2) + 2*k5*covariance15(1,2);
            derivative_k2 = 2*k2*var(error2) + 2*k1*covariance12(1,2) + 2*k3*covariance23(1,2) + 2*k4*covariance24(1,2) + 2*k5*covariance25(1,2);
            derivative_k3 = 2*k3*var(error3) + 2*k2*covariance23(1,2) + 2*k1*covariance13(1,2) + 2*k4*covariance34(1,2) + 2*k5*covariance35(1,2);
            derivative_k4 = 2*k4*var(error4) + 2*k2*covariance24(1,2) + 2*k1*covariance14(1,2) + 2*k3*covariance34(1,2) + 2*k5*covariance45(1,2);
            derivative_k5 = 2*k5*var(error5) + 2*k2*covariance25(1,2) + 2*k1*covariance15(1,2) + 2*k3*covariance35(1,2) + 2*k4*covariance45(1,2);
            
            df = [derivative_k1 derivative_k2 derivative_k3 derivative_k4 derivative_k5];
            
    end


end