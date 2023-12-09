%% Setup
% Define the number of spins and coupling strength
N = 100;  % Number of spins

sigma = sign(randn(N, 1)); % Initialize the vector of spins randomly
sigma(sigma == 0) = 1; % Set zero values to +1

sigma_transpose = sigma'; % Initialize the transpose of sigma

J = (randn(N, N) / sqrt(N));  % Random couplings, properly scaled

% Define the target energy function for the SK model
target_energy = -sum(sum(J .* (sigma' * sigma))) / 2;

% Number of chains and temperatures
num_chains = 19;
temperatures = logspace(-2, 3, num_chains); % Geometric spacing

% Initialize chains
chains = cell(num_chains, 1);
for i = 1:num_chains
    chains{i}.temperature = temperatures(i);
    chains{i}.state = randn(N, 1); % Initialize spins randomly
end

%% Define Functions
% Number of Monte Carlo steps and exchange attempts
num_steps = 2000;
num_exchanges = 400;

% Define the Metropolis-Hastings move function
propose_new_state = @(current_state, temperature) current_state + sqrt(2 * temperature) * randn(N, 1);

% Define compute_energy function
compute_energy = @(chain_state, J) -0.5 * sum(sum(J .* (chain_state' * chain_state)));

% Define accept_state function based on Metropolis criterion
accept_state = @(current_state, proposed_energy, temperature) ...
    rand() < exp(-(proposed_energy - compute_energy(current_state, J)) / temperature);

% Define a function for performing exchange attempts
perform_exchange = @(chain1, chain2) ...
    (rand() < exp((compute_energy(chain1.state, J) - compute_energy(chain2.state, J)) * ...
    (1 / chain1.temperature - 1 / chain2.temperature)));

%% Monte Carlo
% Main PTMC loop
for step = 1:num_steps
    for chain_index = 1:num_chains
        % Perform Metropolis-Hastings move for the current chain
        proposed_state = propose_new_state(chains{chain_index}.state, chains{chain_index}.temperature);

        % Calculate energy for the proposed state
        proposed_energy = compute_energy(proposed_state, J);

        % Accept or reject the proposed state
        if accept_state(chains{chain_index}.state, proposed_energy, chains{chain_index}.temperature)
            chains{chain_index}.state = proposed_state;
        end

        % Ensure spins are +1 or -1
        chains{chain_index}.state = sign(chains{chain_index}.state);
    end

    
   % Perform exchange attempts between chains
    for exchange = 1:num_exchanges
        chain1 = chains{randi(num_chains)};
        chain2 = chains{mod(randi(num_chains) + 1, num_chains) + 1};
        if perform_exchange(chain1, chain2)
            % Swap states between chain1 and chain2
            temp_state = chain1.state;
            chain1.state = chain2.state;
            chain2.state = temp_state;
        end
    end
end

 % Display the final states of each chain
    %disp('Final States:');
    %for i = 1:num_chains
       %disp(['Chain ' num2str(i) ':']);
       %disp(chains{i}.state);
    %end


%% Collect samples from the chains
num_chains = numel(chains);  % Number of chains
num_samples = 100;  % Number of samples to collect from each chain (adjust as needed)

collected_samples = cell(num_chains, num_samples);

% Collect samples from each chain
for i = 1:num_chains
    % Randomly select samples from the stored states of the i-th chain
    selected_indices = randi(size(chains{i}.state, 2), [1, num_samples]);

    % Use a for loop for assignment
    for j = 1:num_samples
        collected_samples{i, j} = chains{i}.state(:, selected_indices(j));
    end
end


%% Data Analysis: selecting samples

% chains{i}.state represents the state/configuration of the i-th chain
num_samples = 100;  % Number of samples for analysis 

% Randomly select samples from each chain for analysis
selected_samples = cell(num_chains, num_samples);
for i = 1:num_chains
    % Randomly select samples from the stored states of the i-th chain
    selected_indices = randi(size(chains{i}.state, 2), [1, num_samples]);
    
    % Use a for loop for assignment
    for j = 1:num_samples
        selected_samples{i, j} = chains{i}.state(:, selected_indices(j));
    end
end

%% Data Analysis: calculate Parisi overlaps, average magnetizations, and spin-spin correlations
N = size(selected_samples{1, 1}, 1);  % Assuming all chains have the same number of spins

parisi_overlaps = zeros(num_chains, num_chains); %initialize Parisi function
for i = 1:num_chains
    for j = 1:num_chains
        if i ~= j
            parisi_overlaps(i, j) = sum(sum(cell2mat(selected_samples(i, :)) .* cell2mat(selected_samples(j, :)))) / (N * num_samples);
        end
    end
end    

% Display the calculated Parisi overlaps
disp('Parisi Overlaps:');
disp(parisi_overlaps);

%num_eigenvalues = 1;  
%largest_magnitude_eigenvalue = eigs(parisi_overlaps, num_eigenvalues, 'largestabs');

% Display the result
%disp('Eigenvalue with Largest Magnitude of Order Parameter Matrix:');
%disp(largest_magnitude_eigenvalue);

% Calculate the eigenvalues of the Parisi overlap matrix
parisi_eigenvalues = eig(parisi_overlaps);

% Display all eigenvalues
disp('Eigenvalues of Parisi Overlap Matrix:');
disp(parisi_eigenvalues);

% Calculate the statistical average of the eigenvalues
average_eigenvalue = mean(parisi_eigenvalues);

% Display the result
disp('Statistical Average of Eigenvalues:');
disp(average_eigenvalue);

% Plot the distribution of eigenvalues
figure;
histogram(parisi_eigenvalues, 'BinEdges', linspace(min(parisi_eigenvalues), max(parisi_eigenvalues), 20), 'Normalization', 'probability');
xlabel('Eigenvalue');
ylabel('Probability Density');
title('Distribution of Eigenvalues of Parisi Overlap Matrix');
grid on;




% Calculate magnetization for each chain
average_magnetizations = zeros(num_chains, 1);

for i = 1:num_chains
    % Calculate average magnetization for each chain
    concatenated_samples = cell2mat(selected_samples(i, :));
    average_magnetizations(i) = mean(concatenated_samples(:));

    % Display the temperature for each chain
    disp(['Chain ' num2str(i) ' Temperature: ' num2str(chains{i}.temperature)]);
end

% Now 'average_magnetizations' is a column vector containing the average magnetization for each chain
disp('Average Magnetizations:');
disp(average_magnetizations);


% Calculate spin-spin correlation for each chain
spin_spin_correlations = zeros(num_chains, 1);

for i = 1:num_chains
    concatenated_samples = cell2mat(selected_samples(i, :));
    %display(transpose(concatenated_samples) .* concatenated_samples)
    spin_spin_correlations(i, 1) = mean(transpose(concatenated_samples) .* concatenated_samples, 'all');
end

% Now 'average_spin_spin_correlations' is a column vector containing the average spin-spin correlation for each chain
disp('Average Spin-Spin Correlations:');
disp(spin_spin_correlations);

%calculate <Si>^2~q
m2=average_magnetizations.*average_magnetizations;
disp('Average Magnetizations^2:');
disp(m2);



%% Figures

% Plot spin-spin correlation (q) on the first figure
figure;
for i = 1:num_chains
    plot(log10(chains{i}.temperature), m2(i), 'g-o', ...
        'DisplayName', ['q Chain ' num2str(i)], 'MarkerSize', 8, 'MarkerFaceColor', 'green', 'LineStyle', '-');
    hold on;
end
ylabel('q');
grid on;
xlabel('Log(Temperature)');
title('Spin-Spin Correlation (q) vs Temperature');
hold off;

% Plot magnetization (m) on the second figure
figure;
for i = 1:num_chains
    plot(log10(chains{i}.temperature), average_magnetizations(i), 'r-o', ...
        'DisplayName', ['Magnetization Chain ' num2str(i)], 'MarkerSize', 8, 'MarkerFaceColor', 'red', 'LineStyle', '-');
    hold on;
end
ylabel('Magnetization');
grid on;
xlabel('Log(Temperature)');
title('Magnetization vs Temperature');
hold off;

% Plot q eig avgs from different runs of this code by log(temp) range 
figure;
qvals= [0 0 0 -3.2761e-17 -2.9121e-17 -6.9161e-17 -1.0192e-16 -3.6401e-17 -2.9121e-17];
Tvals = [5 4.5 1 0.5 0 -0.5 -1 -3 -6];
plot(Tvals, qvals, 'b--o', 'MarkerSize', 8, 'MarkerFaceColor', 'blue');
ylabel('q');
ylim([-1.2E-16 0.2E-16]);
grid on;
xlabel('Log(Temperature)');
title('q vs Temperature');
hold off;
