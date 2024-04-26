% FSU: SENIOR DESIGN TEAM 312
% AUTHORS: THEODOR OWCHARIW, VICTOR BELLERA TOVAR, BEN COVITZ, LUCAS TORRES, LUKE FORBIS
% DATE COMPLETED: 4/25/2024
% Description: Uses a bayesian network and team defined percentages to
%              determine the drowsiness level in a driverthis data 
%              will then be ouptutted via auditory and visua stimuli.

%       Example Bayesian Network Setup

%                 CBR(1)
%               /    |
%              /     |
%             /      |  
%            HP(2)   |     TD(3)
%            |       |       |
%             \      |      /
%              \     |     /
%               \    V    /
%                  D(4) 
% CBR =  breathing rate, true if rate decreases
% HP = Head position, True if head is down
% TD = time duration, True if driving for more than 1 hour 
% D = drowsiness, True if drowsy
addpath('Bayesian Network Framework');
script_dir = fileparts(mfilename('fullpath'));
cd(script_dir);
addpath(genpathKPM(pwd))

N = 5;
if length(ev) > N-1
   error('Number of evidence values exceeds the number of observable nodes')
    return;
end

DUR = 1; CBL = 2; HP = 3; TOD = 4; D = 5;
dag = zeros(N,N);
dag(DUR,CBL) = 1;           % connect DUR to CBL
dag(CBL,[HP D]) = 1;      % connect CBL to HP and D
dag(HP,D) = 1;            % connect HP to D
dag(DUR,D) = 1;             % connect DUR to D
dag(TOD,D) = 1;           % connect TOD to D


node_sizes = [3 2 2 2 2]; % dicrete node sizes 
                          % i.e. CBL and CBH have two values while T and
                          % TOD have 3 possible values


% bnet = mk_bnet(dag, node_sizes, 'observed', onodes);
bnet = mk_bnet(dag, node_sizes);
% CPT tables  false values (1) listed first up to truth values (2) or
% highest int 

CPT_dur = [.6 .25 .15]; % short, med, or long drives
CPT_cbl = [0.3, .2 .15 .7 .8 .85];
CPT_hp = [...
    0.7, 0.3 ... % CBL = 1
    0.8, 0.2, ... % CBL = 2
];
CPT_tod = [.5 .5];  % time of day for driving probs, 7AM-3PM, 3-11 PM, 11PM-7AM
CPT_d = [...
    0.65, 0.55, 0.45, 0.45,	0.4, 0.3, 0.35,	0.3, 0.2, 0.25, 0.15, 0.05, 0.35, ...
    0.4, 0.3, 0.4, 0.35, 0.3, 0.3, 0.2, 0.1, 0.15, 0.08, 0.03, 0.35, 0.45, ...
    0.55, 0.55, 0.6, 0.7, 0.65, 0.7, 0.8, 0.75, 0.85, 0.95, 0.65, 0.6, 0.7, ...
    0.6, 0.65, 0.7,	0.7, 0.8, 0.9, 0.85, 0.92, 0.97];
%% Uncomment for testing purposes
ev(1) = 1; % DUR 1 = short drive, 2 = medium drive, 3 = long drive
ev(2) = 1; % CBL 1 = normal breathing, 2 = slowed breathing
%ev(3) = 1; % HP 1 = head up, 2 = head down
%ev(4) = 1; % TOD 1 = day, 2 = night
bnet.CPD{DUR} = tabular_CPD(bnet, DUR, CPT_dur);  % probability of driving time lengths
bnet.CPD{CBL} = tabular_CPD(bnet, CBL, CPT_cbl);  % probability of breathing rate decreasing 
bnet.CPD{HP} = tabular_CPD(bnet, HP, CPT_hp);     % probability of HP being down, 35% chance the head is down if breathing rate has slowed 
bnet.CPD{TOD} = tabular_CPD(bnet, TOD, CPT_tod);  % probability of TOD during drive
bnet.CPD{D} = tabular_CPD(bnet, D, CPT_d);

% calculate posterior probability 
engine = jtree_inf_engine(bnet);
evidence = cell(1,N);

for i = 1:length(ev)
   evidence{i} = ev(i);
end

[engine, ~] = enter_evidence(engine,evidence);

marg = marginal_nodes(engine,D);
probability_drowsiness=marg.T(2);
probability_drowsiness = probability_drowsiness*100;
disp(probability_drowsiness) % Uncomment for testing
