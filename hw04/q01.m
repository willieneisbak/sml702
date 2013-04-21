function [gSamp,mhSamp] = q01()
    
    % generate data
    mus = [-3,3];
    tau = 10;
    X = normrnd(mus(randi(2,[1,100])),ones(1,100));

    % do Gibbs and MH
    numSamp = 10000;
    gSamp = doGibbs();
    gSamp = gSamp(1000:end,:);
    mhSamp = doMH();
    mhSamp = mhSamp(1000:end,:);

    % ----------
    function gs = doGibbs()
        fprintf('Starting Gibbs\n');
        gs = [];
        mu1 = 0; mu2 = 0;
        z = [ones(1,50),2*ones(1,50)];
        for i=1:numSamp
            X1 = X(find(z==1));
            mu1 = normrnd(sum(X1)/((1/tau)+length(X1)),1/((1/tau)+length(X1)));
            X2 = X(find(z==2));
            mu2 = normrnd(sum(X2)/((1/tau)+length(X2)),1/((1/tau)+length(X2)));
            gs(end+1,:) = [mu1,mu2];
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

    function mhs = doMH()
        fprintf('Starting MH\n');
        mhs = zeros(numSamp,2);
        aRateCount = 0;
        mu = [0,0];
        for i=1:numSamp
            mu_new = mvnrnd(mu,0.1*eye(2));
            num = normpdf(mu_new(1),0,tau)*normpdf(mu_new(2),0,tau)*probXGivenMus(mu_new);
            den = normpdf(mu(1),0,tau)*normpdf(mu(2),0,tau)*probXGivenMus(mu);
            aRatio = num/den; 
            if aRatio>1
                mu = mu_new;
                aRateCount = aRateCount+1;
            else
                u = rand;
                if u < aRatio
                    mu = mu_new; 
                    aRateCount = aRateCount+1;
                end
            end
            mhs(i,:) = mu;
        end
        fprintf('MH acceptance rate = %g\n',aRateCount/numSamp);
    end

    function probX = probXGivenMus(m)
        val1 = 0.5*normpdf(X,m(1),1);
        val2 = 0.5*normpdf(X,m(2),1);
        probX = prod(val1+val2);
    end

end
