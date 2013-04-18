function q01()
    
    % generate data
    mus = [-3,3];
    tau = 10;
    X = normrnd(mus(randi(2,[1,100])),ones(1,100))
    
    % do Gibbs
    gSamp = doGibbs();

    % do MH
    mhSamp = doMH();

    % ----------
    function doGibbs()
        mu1 = 0; mu2 = 0;
        z = [ones(1,50),2*ones(1,50)];
        for i=1:1000
            X1 = X(find(z==1));
            mu1 = normrnd((tau/((1/length(X1))+tau))*mean(X1),1/((1/tau)+n);
            X2 = X(find(z==2));
            mu2 = normrnd((tau/((1/length(X2))+tau))*mean(X2),1/((1/tau)+n);
            for j=1:length(z)
                zPost(1) = normpdf(X(j),mu1,1);
                zPost(2) = normpdf(X(j),mu2,1);
                zPost = zPost/sum(zPost);
                if rand<zPost(1)
                    z(j) = 1;
                else
                    z(j) = 2;
                end
            end
        end
    end

    function doMH()
        %%%%
    end

end
