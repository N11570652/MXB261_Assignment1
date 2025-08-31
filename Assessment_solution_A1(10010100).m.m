%{
(FOR PERSONAL USE TO UNDERSTAND WHAT THE ASSESSMENT WANTS)
MAIN GOAL OF THE ASSESSMENT 
Write MATLAB that simulates N biased random-walkers on a 99×99 grid, 
one at a time, starting on the top row (either fixed at column 50 or 
random along the top). Each walker moves only South/East/West with 
probabilities s, w, e, obeys cyclic E–W boundaries, and stops when 
trying to move South into an occupied cell or the bottom wall. After 
all walkers finish, compute the column heights of the pile and plot a 
histogram. Repeat for the four probability cases and both start-position 
rules, for N = 100 and N = 200. (All directly from the brief.)
%}


clc; clear; close all; 
rng(0); % this is the random seed 



%state the domain size (this sets up the simulation grid) 
rows = 99; 
columns = 99; 


% start setting up the probabilities that have been stated in the task
% sheet 
prob_cases = {[1/3  1/3  1/3], [2/3  1/6  1/6], [4/5  3/10  1/10], [3/5  1/10  3/10]}; %essential creating a vector of s, w and e 


% this part will include the partical counts and the starting positions
N_values = [100, 200];


%the start position 
p_options = {'fixed', 'rand'}; % 'fixed is where all columns start in column 50 and rand is random solumn start (I think....)


% THE OUTER LOOP 
for p_idx = 1:length(p_options) % THis loops over the start positions 
    for n_idx = 1:length(N_values) % tHis part loops over the partical count 
        N = N_values(n_idx);
        P = p_options{n_idx};



        % putting in the new figures for N and P 
        figure; 
        sgtitle(sprintf('random walk: Start=%s, N=%d', P, N));



        for case_idx = 1:length (prob_cases)
            probs = prob_cases{case_idx};
            s = probs(1);
            w = probs(2);
            e = probs(3);


            % run the simulation
            heights = simulate_random_walk(rows, columns, N, P, s, w, e);



            %from the Task sheet the placeholder histogram goes here 
            subplot(2,2, case_idx);
            histogram(heights, 0.5:1:columns+0.5, 'FaceColor', [0.2 0.6 0.8]);
            xlabel('Column index');
            ylabel('height');
            title(sprintf('Case %d, (s =%.2f, w =%.2f, e = %.2f ', case_idx, s, w, e))
        end
    end
end




%simulate the random walk (skim through weeks 1,2,3,4 lecture notes)
function heights = simulate_random_walk(rows, cols, N, P, s, w, e)

lattice = zeros(rows, cols); 

% start looping over the particals. 

for p = 1:N
    if strcmp(P, 'fixed')
        col = 50; % this starts it at the middle 
    else
        col = randi(cols); %if it does not start in the middle it will start at a random column 
    end 

    row = 1; % this indicates the top row 


    while true 
        u = rand;
        if u < s %choosing a rnadom direction based on the probabilties pesented 
            dir = 'S';
        elseif u < s + w
            dir = 'W';
        else 
            dir = 'E';
        end 



    %need to write code that states that if the bottom row is blocked =
    %stop. if the west and east sections are blocked then the simualtion
    %needs to try again.


        if dir == 'S'
            if row == rows || lattice(row+1, col) == 1
                lattice(row, col) = 1;

                break;
            
            else
                row = row + 1;
            end 


        elseif dir == 'W'
            new_col = mod(col-2, cols) + 1;
            if lattice(row, new_col) == 0
                col = new_col;
            end 

        elseif dir == 'E'
            new_col = mod(col, cols) + 1;
            if lattice(row, new_col) == 0
                col = new_col;
            end
        end 
    end
end

heights = sum(lattice,1);
end 

%Part 1 has ended. for purposes of flow structure I will be doing part 2
%within the same sheet. 


%PART 2 



% 1. the true passion distribution 

lambda = 4; 
k = 0:15;
true_pmf = poisspdf(k,lambda);
true_pmf = true_pmf/ sum(true_pmf);  % this will renormalise so it sums to 1 (based on week 3 notes)


%experiment setup 
sample_sizes = [10, 25, 50, 100, 175, 250];
num_experiments = 100;
KL_results = zeros(num_experiments, length(sample_sizes));


%the next thing to do is run the experiment 
for s_idx = 1:length(sample_sizes)
    N = sample_sizes(s_idx);

    for exp = 1:num_experiments
        samples = inverse_transform_poisson(true_pmf, k, N);

        emp_counts = histcounts(samples, -0.5:1:15.5);
        emp_pmf = emp_counts / sum(emp_counts);
        emp_pmf = max(emp_pmf, 1e-10);  %this will avoid the log 


        %the KL deliverence as stated in week 3 notes 
        KL_results(exp, s_idx) = sum(true_pmf .* log(true_pmf ./ emp_pmf));
    end 
end




% mean KL vs the sample size
mean_KL = mean(KL_results,1);
std_KL = std(KL_results, 0, 1);


figure;
errorbar(sample_sizes, mean_KL, std_KL, 'o-', 'linewidth', 1.5);
xlabel('Sample Size');
ylabel('KL Divergence');
title('KL Divergence vs Sample Size');
grid on;


figure;
for s_idx = 1:length(sample_sizes)
    N = sample_sizes(s_idx);
    samples = inverse_transform_poisson(true_pmf, k, N);
    emp_pmf = emp_counts / sum(emp_counts);


    subplot(2,3,s_idx);
    bar(k, [true_pmf; emp_pmf]', 'grouped');
    xlabel('k');
    ylabel('probability');
    title(sprintf('Sample size = %d', N));
    legend('True PMF', 'Emperical PMF');
end







% Inverse sample (from Week 3 lectures)

function samples = inverse_transform_poisson(true_pmf, k, N)

    cdf = cumsum(true_pmf);
    u = rand(1, N);
    samples = zeros(1, N);

for i = 1:N 
    samples(i) = k(find(cdf >= u(i), 1));
end 
end 

















        



    















