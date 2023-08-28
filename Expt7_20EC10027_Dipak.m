%1. Generating a random bit string of size N=256
N = 256;
databits = randi([0,1], 1, N);
%%
%2. Defining a trellis structure for given polynomials.
%Generator polynomials
g1 = [1 1 0];
g2 = [1 1 1];
g3 = [1 0 1];

figure(1)
subplot(211)
stairs(databits);
title("Message Data")

%A trellis should specify the nextstates and output for a given input and current
%state.
numInp = 2;
numOp = 8;
numStates = 4;

nxtState = zeros(numStates, numInp);%It is a numStates x numInp size matrix. 
output = zeros(numStates, numInp);%It is a numStates x numInp size matrix.

for k=1:numStates%current_state
    nxtState(k,1) = floor((k-1)/2);
    nxtState(k,2) = 2 + floor((k-1)/2);
    
    output(k,1) = pow2(2)*xor(xor(g1(1)*0, g1(2)*floor((k-1)/2)),g1(3)*(mod(k-1,2))) + pow2(1)*xor(xor(g2(1)*0, g2(2)*floor((k-1)/2)),g2(3)*(mod(k-1,2))) + pow2(0)*xor(xor(g3(1)*0, g3(2)*floor((k-1)/2)),g3(3)*(mod(k-1,2)));
    output(k,2) = pow2(2)*xor(xor(g1(1)*1, g1(2)*floor((k-1)/2)),g1(3)*(mod(k-1,2))) + pow2(1)*xor(xor(g2(1)*1, g2(2)*floor((k-1)/2)),g2(3)*(mod(k-1,2))) + pow2(0)*xor(xor(g3(1)*1, g3(2)*floor((k-1)/2)),g3(3)*(mod(k-1,2)));
end

%%
%3. Encoding using the trellis structure
state = 0;
codedData = [];

for k=1:N
    temp = de2bi(output(state+1,databits(k)+1), 3);
    codedData = [codedData temp];
    state = nxtState(state+1, databits(k)+1);
end

subplot(212)
stairs(codedData);
title("Encoded data ")

%%
%4. Recieved Data

%CodedData through a Binary symmetric channel
p = 0.3;%probability of error
rcd_data1 = bsc(codedData, p);
%%
%5. Decoding : Vetirbi algorithm
decodedData = viterbiAlgo(rcd_data1, numStates, N, nxtState, output);

%No. of flipped bits and bit error rate.
[count_bits, error_rate] = biterr(decodedData, databits);
%%
%6. Plotting bit error rate as a function of crossover probability.
error = zeros(1, 101);
prob = zeros(1, 101);
for k=0:1000/2
    p_ch = k/1000;
    rcd_data = bsc(codedData, p_ch);
    decoded_data = viterbiAlgo(rcd_data, numStates, N, nxtState, output);
    [bits_flipped, rate] = biterr(decoded_data, databits);
    prob(k+1) = p_ch;
    error(k+1) = rate;
end

figure(2)
plot(prob, error);
title("Bit Error rate Vs crossover probability")
xlim([0 0.6])

%%
% Repeating the same for a binary asymmetric channel

p0 = 0.1;%for symbol '0'
p1 = 0.3;%for symbol '1'

rcd_data2 = basc(codedData, p0, p1);
decodedData2 = viterbiAlgo(rcd_data2, numStates, N, nxtState, output);


%%
function op_data = basc(inp, p0, p1)
    op_data = zeros(1,length(inp));
    for k=1:length(inp)
        if inp(k) == 1
            op_data(k) = bsc(inp(k), p1);
        else
            op_data(k) = bsc(inp(k), p0);
        end
    end
end


function decoded = viterbiAlgo(rcd_data, numStates, N, nxtState, output)
    
    dparr = ones(numStates, N+1)*inf;
    dparr(1,1) = 0;
    %stores the previous element.
    parent = ones(numStates, N+1)*inf;
    parent(1,1) = 0;
    
    for k=1:N
        for j=1:numStates
            if dparr(j,k) ~= inf
                
                dist1 = dparr(j, k) + hammDist(bi2de(rcd_data(3*k-2:3*k)), output(j,1));
                dist2 = dparr(j, k) + hammDist(bi2de(rcd_data(3*k-2:3*k)), output(j,2));
                
                if dist1 < dparr(nxtState(j,1)+1, k+1)
                    dparr(nxtState(j,1)+1, k+1) = dist1;
                    parent(nxtState(j,1)+1, k+1) = j;
                end
                if dist2 < dparr(nxtState(j,2)+1, k+1)
                    dparr(nxtState(j,2)+1, k+1) = dist2;
                    parent(nxtState(j,2)+1, k+1) = j;
                end
                
            end
        end
    end
    
    decoded = [];
    ind = N+1;
    [hamDist,state]  = min(dparr(:, N+1));
    while ind > 1
        decoded = [floor((state-1)/2) decoded];
        state = parent(state,ind);
        ind = ind-1;
    end
end



function hamm_dist = hammDist(num1, num2)
    temp2 = bitxor(num1, num2);
    temp2 = de2bi(temp2);
    hamm_dist = sum(temp2 == 1);
end