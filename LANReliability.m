%% LAN Reliabiltiy
clear;
clc;
close all;

% Defining K, p, and N
K = [1,5,15,50,100];
p = 0:0.01:0.99;
N = 1000;

%% Task 1

% Create sim and calc arrays to be filled
sim = zeros(length(K),length(p));
calc = zeros(length(K),length(p));

% Fill sim and calc arrays
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p)
        pj = p(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runSingleLinkSim(Ki,pj,N);

        % Use calculation and fill calculation array
        calc(i,j) = Ki/(1 - pj);
    end
end

% Plot results
for i = 1:length(K)
    Ki = K(i);
    
    figure;
    semilogy(p,calc(i,:),'color','r')
    hold on;
    semilogy(p,sim(i,:),'ko','color','b');
    title(sprintf('Task 1 K = %d', Ki));
    xlabel('Probability');
    ylabel('Number of Transmissions');
    legend('Calculated','Simulation');
    hold off;
end

% Plot results on one graph
figure;
semilogy(p,sim(1,:),'ko','color','b');
hold on;
semilogy(p,sim(2,:),'ko','color','r');
semilogy(p,sim(3,:),'ko','color','g');
semilogy(p,sim(4,:),'ko','color','c');
semilogy(p,sim(5,:),'ko','color','m');
semilogy(p,calc(1,:),'color','k');
semilogy(p,calc(2,:),'color','k');
semilogy(p,calc(3,:),'color','k');
semilogy(p,calc(4,:),'color','k');
semilogy(p,calc(5,:),'color','k');
legend('K=1','K=5','K=15','K=50','K=100');
title('Task 1 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% Task 2

% Create sim and calc arrays to be filled
sim = zeros(length(K),length(p));
calc = zeros(length(K),length(p));

% Fill sim and calc arrays
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p)
        pj = p(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runTwoSeriesLinkSim(Ki,pj,N);

        % Use calculation and fill calculation array
        calc(i,j) = Ki/(1 - pj)^2;
    end
end

% Plot results
for i = 1:length(K)
    Ki = K(i);
    
    figure;
    semilogy(p,calc(i,:),'color','r')
    hold on;
    semilogy(p,sim(i,:),'ko','color','b');
    title(sprintf('Task 2 K = %d', Ki));
    xlabel('Probability');
    ylabel('Number of Transmissions');
    legend('Calculated','Simulation');
    hold off;
end

% Plot results on one graph
figure;
semilogy(p,sim(1,:),'ko','color','b');
hold on;
semilogy(p,sim(2,:),'ko','color','r');
semilogy(p,sim(3,:),'ko','color','g');
semilogy(p,sim(4,:),'ko','color','c');
semilogy(p,sim(5,:),'ko','color','m');
semilogy(p,calc(1,:),'color','k');
semilogy(p,calc(2,:),'color','k');
semilogy(p,calc(3,:),'color','k');
semilogy(p,calc(4,:),'color','k');
semilogy(p,calc(5,:),'color','k');
legend('K=1','K=5','K=15','K=50','K=100');
title('Task 2 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% Task 3

% Create sim array to be filled
sim = zeros(length(K),length(p));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p)
        pj = p(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runTwoParallelLinkSim(Ki,pj,N);
    end
end

% Plot results
for i = 1:length(K)
    Ki = K(i);
    
    figure;
    semilogy(p,sim(i,:),'ko','color','b');
    title(sprintf('Task 3 K = %d', Ki));
    xlabel('Probability');
    ylabel('Number of Transmissions');
end

% Plot results on one graph
figure;
semilogy(p,sim(1,:),'ko','color','b');
hold on;
semilogy(p,sim(2,:),'ko','color','r');
semilogy(p,sim(3,:),'ko','color','g');
semilogy(p,sim(4,:),'ko','color','c');
semilogy(p,sim(5,:),'ko','color','m');
legend('K=1','K=5','K=15','K=50','K=100');
title('Task 3 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% Task 4

% Create sim array to be filled
sim = zeros(length(K),length(p));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p)
        pj = p(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCompoundNetworkSim(Ki,pj,N);
    end
end

% Plot results
for i = 1:length(K)
    Ki = K(i);
    
    figure;
    semilogy(p,sim(i,:),'ko','color','b');
    title(sprintf('Task 4 K = %d', Ki));
    xlabel('Probability');
    ylabel('Number of Transmissions');
end

% Plot results on one graph
figure;
semilogy(p,sim(1,:),'ko','color','b');
hold on;
semilogy(p,sim(2,:),'ko','color','r');
semilogy(p,sim(3,:),'ko','color','g');
semilogy(p,sim(4,:),'ko','color','c');
semilogy(p,sim(5,:),'ko','color','m');
legend('K=1','K=5','K=15','K=50','K=100');
title('Task 4 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% Task 5

% Initialize K values
K = [1,5,10];

%% 5.1

% Initialize probabilities
p1 = 0.1;
p2 = 0.6;
p3 = 0:0.01:0.99;

% Create sim array to be filled
sim = zeros(length(K),length(p3));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p3)
        p3j = p3(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1,p2,p3j,N);
    end
end

% Plot results on one graph
figure;
semilogy(p3,sim(1,:),'ko','color','b');
hold on;
semilogy(p3,sim(2,:),'ko','color','r');
semilogy(p3,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.1 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% 5.2

% Initialize probabilities
p1 = 0.6;
p2 = 0.1;
p3 = 0:0.01:0.99;

% Create sim array to be filled
sim = zeros(length(K),length(p3));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p3)
        p3j = p3(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1,p2,p3j,N);
    end
end

% Plot results on one graph
figure;
semilogy(p3,sim(1,:),'ko','color','b');
hold on;
semilogy(p3,sim(2,:),'ko','color','r');
semilogy(p3,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.2 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% 5.3

% Initialize probabilities
p1 = 0.1;
p2 = 0:0.01:0.99;
p3 = 0.6;

% Create sim array to be filled
sim = zeros(length(K),length(p2));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p2)
        p2j = p2(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1,p2j,p3,N);
    end
end

% Plot results on one graph
figure;
semilogy(p2,sim(1,:),'ko','color','b');
hold on;
semilogy(p2,sim(2,:),'ko','color','r');
semilogy(p2,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.3 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% 5.4

% Initialize probabilities
p1 = 0.6;
p2 = 0:0.01:0.99;
p3 = 0.1;

% Create sim array to be filled
sim = zeros(length(K),length(p2));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p2)
        p2j = p2(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1,p2j,p3,N);
    end
end

% Plot results on one graph
figure;
semilogy(p2,sim(1,:),'ko','color','b');
hold on;
semilogy(p2,sim(2,:),'ko','color','r');
semilogy(p2,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.4 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% 5.5

% Initialize probabilities
p1 = 0:0.01:0.99;
p2 = 0.1;
p3 = 0.6;

% Create sim array to be filled
sim = zeros(length(K),length(p1));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p1)
        p1j = p1(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1j,p2,p3,N);
    end
end

% Plot results on one graph
figure;
semilogy(p1,sim(1,:),'ko','color','b');
hold on;
semilogy(p1,sim(2,:),'ko','color','r');
semilogy(p1,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.5 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;

%% 5.6

% Initialize probabilities
p1 = 0:0.01:0.99;
p2 = 0.6;
p3 = 0.1;

% Create sim array to be filled
sim = zeros(length(K),length(p1));

% Fill sim array
for i = 1:length(K)
    Ki = K(i);
    for j = 1:length(p1)
        p1j = p1(j);
        
        % Run simulation and fill simulation array
        sim(i,j) = runCustomCompoundNetworkSim(Ki,p1j,p2,p3,N);
    end
end

% Plot results on one graph
figure;
semilogy(p1,sim(1,:),'ko','color','b');
hold on;
semilogy(p1,sim(2,:),'ko','color','r');
semilogy(p1,sim(3,:),'ko','color','g');
legend('K=1','K=5','K=10');
title('Task 5.6 All K Values');
xlabel('Probability');
ylabel('Number of Transmissions');
hold off;