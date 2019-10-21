function simNetworks()
% simNetworks
% Simulates and compares CEBCRA, LEACH routing protocols for a given
% wireless sensor network. 
% 
% Functions:
% initializeNetwork
% setupCEBCRA
%
% Change Log
%  Function definitions, node parameters
%  Added plotNetworks
%  Updated plotNetworks


% parameters setup - these could be encapsulated in node function
% n is network structure for common parameters
n.numNodes=150; %500
n.numActiveNodes=n.numNodes;
n.numDeadNodes=0;
n.rxRate=10^6; %RX rate 1 Mbps
n.messageSize=500; %Byte, packet message size
n.packetHeaderSize=25; %Byte,packet header size
n.initialEnergy=2; %Joule, initial energy per node
n.maxDist=25; %25,30 neighbor dist, calc as max distance at fixed transmission energy
n.colorsCl=rand(180,3);
n.packetsToBS(1)=0;
n.activeNodes(1)=n.numNodes;
n.consumedEnergy(1)=0;
n.roundToNewCH=20;


global MAXDIST;
MAXDIST=n.maxDist;

%initialize networks, nodes,
%same network used for all three 
n=initializeNetwork(n); %node locations, parameters, base station
nCEBCRA=n;

nLEACH=n;
nLEACH=initializeNetworkLEACH(nLEACH);
clear n;

%separate figure number for each network
nCEBCRA.figNum=100; 

nLEACH.figNum=302;


%select which protocol to run
%uncomment below, 1 at a time

round=1;
while 1
    disp(round)

    
    %CEBCRA
    if (round==1)||mod(round,nCEBCRA.roundToNewCH)==0
    nCEBCRA=setupCEBCRA(nCEBCRA);
    end
    nCEBCRA=steadyStateCEBCRA(nCEBCRA);
    plotNetwork(nCEBCRA);
    if (nCEBCRA.numDeadNodes>nCEBCRA.numNodes-10)
        break; 
    end
    


     %LEACH
     if (round==1)||mod(round,nLEACH.roundToNewCH)==0
         nLEACH=setupLEACH(nLEACH,round);
     end
     nLEACH=steadyStateLEACH(nLEACH,round);
     plotNetwork(nLEACH);
     if (nLEACH.numDeadNodes>nLEACH.numNodes-10)
         break; 
     end
    


    
    round=round+1;
end
    
    
% calculate and saved results
if exist('nLEACH500.mat','file')
load('nLEACH500.mat');
x=length(fieldnames(n));
n.(['n' num2str(x+1)])=nLEACH;
save('nLEACH500.mat','n');
else
    n.n1=nLEACH;
save('nLEACH500.mat','n');
end
%calculateResults();
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=initializeNetwork(n)
%initalize networks
x_max=200; %meters, square field
y_max=200; %meters, square field


%initialize nodes
for i=1:n.numNodes
    
    %location
    n.node(i).x=x_max*rand;
    n.node(i).y=y_max*rand;
    
    %initial energy
    n.node(i).energy=n.initialEnergy;
    
    % node type
    n.node(i).type='n'; %n- node, CH-cluster head
    
    % node status
    n.node(i).status='A'; %A-Alive or D-Dead
    
    % sda
end

%initialize base station
n.BS.x=175; %275 meters for 500
n.BS.y=225; %325meters


function n=initializeNetworkLEACH(n)
%extra parameters for LEACHinitalize networks

x_max=200; %meters, square field
y_max=200; %meters, square field
%initialize nodes
for i=1:n.numNodes

    n.node(i).C_i=1;
    n.node(i).C_i_counter=0;
    n.node(i).d_toCH=3*x_max;
    n.node(i).numNeighbors=0;

end



function plotNetwork(n)
%plots alive and dead Nodes, BS,
% add energy and data packets sent plot
%possibly animate message sent for demo

colorsCl=n.colorsCl;
% colorsCl(1,:)=[0 0 0];
% for i=2:30
%     colorsCl(i,:)=colorsCl(i-1,:)+[1/35 1/35 0];
% end

figure(n.figNum); hold off;
plot(1,1); hold on;
xlim([-30 n.BS.x+30]);
ylim([-30 n.BS.y+30]);

%Cluster heads to base station plot
plot(n.BS.x,n.BS.y,'gX');
plot(n.BS.x,n.BS.y,'gO');
for i=1:n.numNodes
    iCH=find(n.CH==n.node(i).CH);
    ncolor=colorsCl(iCH,:);
    %If Node is alive, color is blue; if dead, red
    if strcmp(n.node(i).status,'A')
        % nodes are 'x's and CH are 'o's
        if strcmp(n.node(i).type,'CH')
            plot(n.node(i).x,n.node(i).y,'bo','linewidth',2);%,'color',ncolor);
            plot(n.node(i).x,n.node(i).y,'bx','linewidth',2);%,'color',ncolor);
        else
            plot(n.node(i).x,n.node(i).y,'bo','color',ncolor);
        end
    else
        plot(n.node(i).x,n.node(i).y,'rx','color',[1 0 0]);
    end
end
hold off;
% xlabel('X (m)','FontSize',12,'FontWeight','bold')
% ylabel('Y (m)','FontSize',12,'FontWeight','bold')

figure(n.figNum+1); 
subplot(2,1,1); hold on;
plot(n.packetsToBS)
xlabel('Round')
ylabel('Packets to BS')
subplot(2,1,2); hold on;
plot(n.activeNodes)
xlabel('Round')
ylabel('Active Nodes')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CEBRCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=setupCEBCRA(n)
%Setup cluster head selection using TDMA, cluster setup 
%input: Residual Energy nodes 
%output: Nodes as CH or Cluster member

% clear previous clusterheads
n.CH=[];

%setp 1: node i transmits HELLO message to all sensor nodes in range
%also clear old costs, weights
for i=1:n.numNodes
    n.node(i).weight=0;
    n.node(i).cost=[];
    n.node(i).neighbor=[];
    n.node(i).CH=0;
    
    % if node doesn't have enough energy to TX, node is dead
    if (n.node(i).energy<.09)
        
        if strcmp(n.node(i).status,'A')
            n.node(i).status='D';
            n.numDeadNodes=n.numDeadNodes+1;
        end
        
    else
        
        %step 2: node i counts all received HELLO messages and calcs W(n)
        nCount=0;
        for ii=1:n.numNodes
            if (i~=ii)&&(sqrt((n.node(i).x-n.node(ii).x)^2+(n.node(i).y-n.node(ii).y)^2)<n.maxDist)
                nCount=nCount+1;
                n.node(i).neighbor(nCount)=ii;
            end
        end
        n.node(i).energy=n.node(i).energy-transmitEnergy(); %fixed power Tx
        n.node(i).energy=n.node(i).energy-nCount*receiveEnergy();
        
        
        n.node(i).weight=(n.node(i).energy)^2/nCount;
    end
end

%step 3: node i broadcast W(n) and receives W(n) for all neighbors
%step 4: set flag
%step 5: Receive neighbors weights and compare to W(n), if neighbor > node
%i gives up CH
%step 6: if flag==1, declare self as CH, recieve join requests, else
%recieve join advertise rquest from CH, cal cost of each CH and CH with
%highest cost
nCH=1;
for i=1:n.numNodes
    if strcmp(n.node(i).status,'A')
        flag=1;
        for ii=1:length(n.node(i).neighbor)
            iNeighbor=n.node(i).neighbor(ii);
            if (n.node(i).weight<n.node(iNeighbor).weight)
                n.node(i).type='n'; flag=0; break;
            end
        end
        %declare CH and
        if flag
            n.node(i).type='CH';
            n.node(i).CH=i;
            n.node(i).msgBuffer=0; %num of message in Buffer
            n.CH(nCH)=i;
            nCH=nCH+1;
        end
    end
end

for i=1:n.numNodes
    if strcmp(n.node(i).status,'A')
    if strcmp(n.node(i).type,'CH')
        %disp('Broadcast join request');
        n.node(i).energy=n.node(i).energy-transmitEnergy(); %fixed power Tx
    else
        distBS=sqrt((n.node(i).x-n.BS.x)^2+(n.node(i).y-n.BS.y)^2);
        for ii=1:length(n.CH)
            iiCH=n.CH(ii);
            distCH=sqrt((n.node(i).x-n.node(iiCH).x)^2+(n.node(i).y-n.node(iiCH).y)^2);
            n.node(i).cost(ii)=n.node(i).energy*distBS/distCH;
        end
        [cMax iCH]=max(n.node(i).cost);
        n.node(i).CH=n.CH(iCH);
    end
    end
end
if length(n.CH)<3
x=1;
end


function n=steadyStateCEBCRA(n)
%steady state operation for CEBCRA
%input: residual energies of all CH from neighbors
%       distance C to BS and all Ch from neightbors
%output: next hop relay node of C

% TDMA Schedule section - simplified for simulation
% assume each cluster is well behaved and no error in TX/RX
%receive data from local cluster
for i=1:n.numNodes
    if strcmp(n.node(i).status,'A')
    %send message to CH
    if ~strcmp(n.node(i).type,'CH')
        iCH=n.node(i).CH; %CH index of node i
        n.node(i).energy=n.node(i).energy-transmitEnergy(n.node(i),n.node(iCH));
        %n.node(i).energy=n.node(i).energy-transmitEnergy();
        n.node(iCH).energy=n.node(iCH).energy-receiveEnergy();
        n.node(iCH).msgBuffer=n.node(iCH).msgBuffer+1;
    else
        n.node(i).msgBuffer=n.node(i).msgBuffer+1;
    end
    
    end
end

%aggregate data
for i=1:length(n.CH)
    iCH=n.CH(i);
    n.node(iCH).energy=n.node(iCH).energy-aggEnergy(n.node(iCH).msgBuffer);
end

% select relay CH node 
minDist=10000; closestCH=1;%n.CH(1);
for i=1:length(n.CH)
    relay=0;%node 0 is BS
    iCH=n.CH(i);
    %n.node(iCH).relay=re;%node 0 is BS
    costRelay=1/transmitEnergy(n.node(iCH),n.BS);
    dBS=distance(n.node(iCH),n.BS);
    if dBS<minDist; minDist=dBS; closestCH=iCH;end
    if dBS<2*n.maxDist
    for ii=1:length(n.CH)
        if i~=ii
            iiCH=n.CH(ii);
            costRelayii=n.node(iiCH).energy/(transmitEnergy(n.node(iCH),n.node(iiCH))* ...
                (transmitEnergy(n.node(iiCH),n.BS)+receiveEnergy()));
                %dRtoBS=distance(n.node(iiCH),n.BS); 
            if (costRelayii>costRelay)%&&(dRtoBS<dBS)
                relay=iiCH;
                costRelay=costRelayii;
                %if self is selected, tx to BS
                if i==ii
                    relay=0;
                end
            end
        end
    end
    end

    n.node(iCH).relay=relay;
end
n.node(closestCH).relay=0;

%relay data to BS
relayToBS=1;
packetsToBS=0;
rcount=0;
while relayToBS
    %for each CH, relay msg to next hop
    numTX=0;
    for i=1:length(n.CH)
        iCH=n.CH(i);
        numMsg=n.node(iCH).msgBuffer;
        if numMsg~=0
        %TX messages to next hop relay
        %for ii=1:length(n.node(iCH).msgBuffer)
            if n.node(iCH).relay==0 
                % BS is relay 0, send msg directly to BS
                n.node(iCH).energy=n.node(iCH).energy-transmitEnergy(n.node(iCH),n.BS);%add numMsg*
                %n.node(iCH).energy=n.node(iCH).energy-numMsg*transmitEnergyCEBCRA(n.node(iCH),n.BS);
                n.node(iCH).msgBuffer=0;
                packetsToBS=packetsToBS+numMsg;
            else
                % relay message to relay node,
                iRelay=n.node(iCH).relay;
                %n.node(iCH).energy=n.node(iCH).energy-numMsg*transmitEnergyCEBCRA(n.node(iCH),n.node(iRelay));
                n.node(iCH).energy=n.node(iCH).energy-transmitEnergy(n.node(iCH),n.node(iRelay));%add numMsg*
                n.node(iRelay).energy=n.node(iRelay).energy-numMsg*receiveEnergy();
                n.node(iCH).msgBuffer=0;
                n.node(iRelay).msgBuffer=n.node(iRelay).msgBuffer+numMsg;
            end
        else
            numTX=numTX+1;
        end
            
       % end
        
    end
    rcount=rcount+1;
    % if all TX received, end steady state
    if numTX==length(n.CH); break; end
    if rcount>1000; break; end
    
end

%update protocol performance measurements
n.packetsToBS(length(n.packetsToBS)+1)=n.packetsToBS(length(n.packetsToBS))+packetsToBS;
n.activeNodes(length(n.activeNodes)+1)=n.numNodes-n.numDeadNodes;
sumE=0;
for i=1:n.numNodes
   if (n.node(i).energy)>0
    sumE=sumE+n.node(i).energy;
    end
end
n.consumedEnergy(length(n.consumedEnergy)+1)=(n.initialEnergy*n.numNodes-sumE);
%n.consumedEnergy(length(n.consumedEnergy)+1)=n.consumedEnergy(length(n.consumedEnergy))+(n.initialEnergy*n.numNodes-sumE);
%n.consumedEnergy(length(n.consumedEnergy))+(n.initialEnergy*n.numNodes-sumE)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LEACH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n=setupLEACH(n,r)
N=n.activeNodes(r);
if N~=0
M=200;
epsilon_fs=10*10^-12; %J/bit/m^2
epsilon_mp=0.0013*10^-12; %J/bit/m^4
d_toBS=sqrt((175-M/2)^2+(225-M/2)^2);
k_opt=(N*epsilon_fs*M^2/(2*pi*epsilon_mp*d_toBS^4))^0.5;
k=k_opt;
%{
% In this paper,
%we assume these parameters are programmed into the nodes a
%priori. However, this approach does not work well in dynamic
%networks.

E_CH=l*E_elec*(N/k-1)+l*E_DA*N/k+l*E_elec+l*epsilon_mp*d_toBS^2
E_nonCH=l*E_elec+l*epsilon_fs*d_toCH^2
E_nonCH=l*E_elec+l*epsilon_fs*(M^2/(2*pi*k))^2

k_opt=(N*epsilon_fs*M^2/(2*pi*epsilon_mp*d_toBS^4))^0.5
%}
E_total=0;
for i=1:n.numNodes
    if strcmp(n.node(i).status,'A')
        E_total=E_total+n.node(i).energy;
    end
end
for i=1:n.numNodes
    if (n.node(i).energy<.2)&&strcmp(n.node(i).status,'A')
        n.node(i).status='D'; n.numDeadNodes=n.numDeadNodes+1;
    end
    if strcmp(n.node(i).status,'A')
        % Assign to nodes probability of becoming a CH.
        n.node(i).P_i=min(k*n.node(i).energy/E_total,1);
        %{
        if n.node(i).C_i==1
            n.node(i).P_i=1/(N/k-mod(r,N/k)); % r is round number
        elseif n.node(i).C_i==0
            n.node(i).P_i=0;
        end
        %}
    end
end
% Select CHs
nCH=0;
n.CH=[];
w=1;
while nCH==0 %ensure at least one CH exists.
    count=0;
    for i=1:n.numNodes
        n.node(i).CH=0;
        if strcmp(n.node(i).status,'A')
            count=count+1;
            if rand<n.node(i).P_i
                n.node(i).type='CH';
                n.node(i).CH=i;
                nCH=nCH+1;
                n.CH(nCH)=i;
                n.node(i).C_i=0;
                n.node(i).C_i_counter=mod(r,N/k); 
                n.node(i).energy=n.node(i).energy-transmitEnergyLEACH(); %advertisement message via non-persistent CSMA
            else
                n.node(i).type='n';
            end
        end
    end
    if count==0 %no alive nodes
        break
    end
    w=w+1;
    if w>10000
        %no nodes with enough energy to be CH for next round, return
        n.numDeadNodes=n.numNodes;
        xxxx=1;
        return;
    end
end
% Update counters for determining probabilities of becoming a CH.
for i=1:n.numNodes
    n.node(i).C_i_counter=max(n.node(i).C_i_counter-1,0);
    if n.node(i).C_i_counter==0
        n.node(i).C_i=1;
    end
end
try
% Nodes choose CHs
for i=1:n.numNodes %for each node
    n.node(i).d_toCH=300;
    if strcmp(n.node(i).status,'A') %that's alive
        for ii=1:n.numNodes % iterate through
            %ii
            if (strcmp(n.node(ii).status,'A'))&&(strcmp(n.node(ii).type,'CH'))&&(ii~=i) %alive CHs excluding self,
                d_toCH2=sqrt((n.node(i).x-n.node(ii).x)^2+(n.node(i).y-n.node(ii).y)^2); %calculate distance to them,
                n.node(i).energy=n.node(i).energy-receiveEnergyLEACH(1); %(receives signal from all CHs)
                %n.node(i).d_toCH
                if d_toCH2<n.node(i).d_toCH %if currently its the closest, %have to initialize to 2*M?
                    n.node(i).d_toCH=d_toCH2; %update minimum distance
                    n.node(i).CH=ii; %and select as new CH.
                end
            end
        end
        n.node(n.node(i).CH).numNeighbors=n.node(n.node(i).CH).numNeighbors+1; %increment size of CH's neighbors list %have to initialize at 0?
        n.node(n.node(i).CH).neighbors(n.node(n.node(i).CH).numNeighbors)=i; %append node ID to CH's neighbors list
        n.node(i).energy=n.node(i).energy-transmitEnergyLEACH(n.node(i),n.node(n.node(i).CH),1); %decrease energy from transmission from node to CH
        n.node(n.node(i).CH).energy=n.node(n.node(i).CH).energy-receiveEnergyLEACH(1); %decrease energy from receiving
    end
end
catch
    dd=2;
end
end

function n=steadyStateLEACH(n,r)
if n.numNodes==n.numDeadNodes;return;end
N=n.activeNodes(r);
if N~=0
for i=1:n.numNodes %for each node
    if strcmp(n.node(i).status,'A')&&(strcmp(n.node(i).type,'CH')) %that's an alive CH
        n.node(i).energy=n.node(i).energy-transmitEnergyLEACH(); %set up TDMA schedule and transmit it
        for ii=1:n.node(i).numNeighbors
            n.node(n.node(i).neighbors(ii)).energy=n.node(n.node(i).neighbors(ii)).energy-receiveEnergyLEACH(1); %receive TDMA schedule
            n.node(n.node(i).neighbors(ii)).energy=n.node(n.node(i).neighbors(ii)).energy-transmitEnergyLEACH(n.node(i),n.node(n.node(i).neighbors(ii)),0); %transmit sensor data
        end
        n.node(i).energy=n.node(i).energy-receiveEnergyLEACH(0)*n.node(i).numNeighbors; %receives data from sensors
        n.node(i).energy=n.node(i).energy-aggregateDataEnergyLEACH(n.node(i).numNeighbors); %aggregate data
        n.node(i).energy=n.node(i).energy-transmitEnergyLEACH(n.node(i),n.BS,0); %and send to BS
    end
end
packetsToBS=N;
n.packetsToBS(length(n.packetsToBS)+1)=n.packetsToBS(length(n.packetsToBS))+packetsToBS;
n.activeNodes(length(n.activeNodes)+1)=n.numNodes-n.numDeadNodes;
sumE=0;
for i=1:n.numNodes
    if (n.node(i).energy)>0
    sumE=sumE+n.node(i).energy;
    end
end
n.consumedEnergy(length(n.consumedEnergy)+1)=(n.initialEnergy*n.numNodes-sumE);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%DEBR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function n=setupDEBR(n)
%step 1: node i determines power to base station

for i=1:n.numNodes
    n.node(i).ECik=[];      %Energy Cost node i to neighbor node k
    n.node(i).ECkbs=[];     %Energy Cost node k to BS
    n.node(i).Kenergy=[];   %Energy at node k
    n.node(i).TEC=[];       %Total Energy Cost for each path node i to BS
    n.node(i).neighbor=[];
    n.node(i).distBS=distance(n.node(i),n.BS);
    n.node(i).ECbs=transmitEnergy(n.node(i),n.BS,500);%energy from node i to base station
    
    % if node doesn't have enough energy to TX, node is dead
    if (n.node(i).energy<.01)
        if strcmp(n.node(i).status,'A')
            n.node(i).status='D';
            n.numDeadNodes=n.numDeadNodes+1;
        end
        
    else
        %step 2: node i counts all received SETUP messages and distance to
        %neighbors
        nCount=0;
        for ii=1:n.numNodes
            if (i~=ii)&&(sqrt((n.node(i).x-n.node(ii).x)^2+(n.node(i).y-n.node(ii).y)^2)<n.maxDist)
                nCount=nCount+1;
                n.node(i).neighbor(nCount)=ii;
                n.node(i).ECik(nCount)=transmitEnergy(n.node(i),n.node(ii),500);
                %ECik is Energy Cost from node i to neighbor node k
            end
        end
        n.node(i).energy=n.node(i).energy-transmitEnergy(); %fixed power Tx
        n.node(i).energy=n.node(i).energy-nCount*receiveEnergy();
    end
    
end
%step 3: Initialize Energy Cost Tables
for i=1:n.numNodes
    n.node(i).TEC(1)=n.node(i).ECbs/n.node(i).energy;
    for ii=1:size(n.node(i).neighbor,2)
        n.node(i).ECkbs(ii)=n.node(n.node(i).neighbor(ii)).ECbs;
        n.node(i).Kenergy(ii)=n.node(n.node(i).neighbor(ii)).energy;
        n.node(i).TEC(ii+1)=(n.node(i).ECik(ii)/n.node(i).energy)+(n.node(i).ECkbs(ii)/n.node(i).Kenergy(ii));
    end
end

%function n=steadyStateDEBR(n)
%steady state operation for DEBR
%input: residual energies of all neighbors
%Node i selects node k as best candidate to send data
%If node i is best candidate, or ties for best candidate, send data direct
%to BS.
%If multiple best candidates exist, arbitrarily pick one
%update TEC tables after each Tx.
%A=[];
%for i=1:n.numNodes
%    [~,I]=min(n.node(i).TEC);
    
%    if I~=1
%        A(1,i)=n.node(i).neighbor(I-1);
%        [~,II]=min(n.node(n.node(i).neighbor(I-1)).TEC);
%        if II~=1
%            A(2,i)=n.node(n.node(i).neighbor(I-1)).neighbor(II-1);
%        else
%            A(2,i)=0;
%        end    
%    else    
%        A(1,i)=0;
%        A(2,i)=0;    
%    end   
%    if i-A(2,i)==0
%        A(3,i)=888;
%    end    
%end

loopCount=0; %debug
t2flag=0;t3flag=0;%debug
for i=1:n.numNodes
    % if node doesn't have enough energy to TX, node is dead
    if (n.node(i).energy<.25)
        if strcmp(n.node(i).status,'A')
            n.node(i).status='D';
            n.node(i).ECik=ones(size(n.node(i).ECik));
            n.numDeadNodes=n.numDeadNodes+1;
        end
        
    else
    if strcmp(n.node(i).status,'A')
        [~,I]=min(n.node(i).TEC);   %find position of lowest TEC for node i
        if I~=1                     %If Tx to node k, subtract energy from node i to node k with lowest TEC
            t=i;t1=0;t2=0;%debug
                                    %  [~,I]=min(n.node(t).TEC);
            while I~=1
               try
                loopCount=loopCount+1;
                n.node(t).energy=n.node(t).energy-n.node(t).ECik(I-1); %reduce node i by ECik
                %For each neighbor node k, reduce energy for Rx
                %For node i, update new Kenergy & TEC values
                for ii=1:size(n.node(t).neighbor,2) 
                    n.node(n.node(t).neighbor(ii)).energy=n.node(n.node(t).neighbor(ii)).energy-receiveEnergy();
                    n.node(t).Kenergy(ii)=n.node(n.node(t).neighbor(ii)).energy;
                    n.node(t).TEC(ii+1)=(n.node(t).ECik(ii)/n.node(t).energy)+(n.node(t).ECkbs(ii)/n.node(t).Kenergy(ii));
                    %For each neighbor k, update Kenergy & TEC values
                    for iii=1:size(n.node(n.node(t).neighbor(ii)).neighbor,2)
                        n.node(n.node(t).neighbor(ii)).Kenergy(iii)=n.node(t).energy;
                        n.node(n.node(t).neighbor(ii)).TEC(iii+1)= ...
                            (n.node(n.node(t).neighbor(ii)).ECik(iii)/n.node(n.node(t).neighbor(ii)).energy)+...
                            (n.node(n.node(t).neighbor(ii)).ECkbs(iii)/n.node(n.node(t).neighbor(ii)).Kenergy(iii));
                    end 
                end
    %            t2=t1;
    %            t1=t;
                t=n.node(t).neighbor(I-1);
                [~,I]=min(n.node(t).TEC);
    %            if t2==t
    %                t2flag=t2flag+1;t2a=t2;
    %                if n.node(t).TEC(1)-min(n.node(t).TEC)<0.0001
    %                    I=1;t3flag=t3flag+1;
    %                else   
    %                    [~,I]=min(n.node(t).TEC);
    %                end
    %            end
               catch errorMsg
                   disp(errorMsg)
                   n.node(i).energy=n.node(i).energy-n.node(i).ECbs;
                   break;
               end
               if loopCount>500;
                   n.node(i).energy=n.node(i).energy-n.node(i).ECbs;
                   break;
                   
               end
           end    
       else      %If lowest TEC is direct to BS,subtract energy to BS
           try
           n.node(i).energy=n.node(i).energy-n.node(i).ECbs;
           %For each neighbor node k, reduce energy for Rx
                %For node i, update new Kenergy & TEC values
                for ii=1:size(n.node(i).neighbor,2) 
                    n.node(n.node(i).neighbor(ii)).energy=n.node(n.node(i).neighbor(ii)).energy-receiveEnergy();
                    n.node(i).Kenergy(ii)=n.node(n.node(i).neighbor(ii)).energy;
                    n.node(i).TEC(ii+1)=(n.node(i).ECik(ii)/n.node(i).energy)+(n.node(i).ECkbs(ii)/n.node(i).Kenergy(ii));
                    %For each neighbor k, update Kenergy & TEC values
                    for iii=1:size(n.node(n.node(i).neighbor(ii)),2)
                        n.node(n.node(i).neighbor(ii)).Kenergy(iii)=n.node(i).energy;
                        n.node(n.node(i).neighbor(ii)).TEC(iii+1)= ...
                            (n.node(n.node(i).neighbor(ii)).ECik(iii)/n.node(n.node(i).neighbor(ii)).energy)+...
                            (n.node(n.node(i).neighbor(ii)).ECkbs(iii)/n.node(n.node(i).neighbor(ii)).Kenergy(iii));
                    end 
                end
      %          for ii=1:size(n.node(i).neighbor,2) %Subtract Rx energy for each neighbor k
      %          n.node(n.node(i).neighbor(ii)).energy=n.node(n.node(i).neighbor(ii)).energy-receiveEnergy();
      %          %For each neighbor k, update Kenergy & TEC for Tx node i
      %          for iii=1:size(n.node(n.node(i).neighbor(ii)),2)
      %              if n.node(n.node(i).neighbor(ii)).neighbor(iii)==i
      %                 n.node(n.node(i).neighbor(ii)).Kenergy(iii)=n.node(i).energy;
      %                 n.node(n.node(i).neighbor(ii)).TEC(iii+1)=...
      %                 (n.node(i).ECik(iii)/n.node(i).energy)+...
      %                 (n.node(i).ECkbs(iii)/n.node(i).Kenergy(iii));
      %              end
      %          end
       %    end
           catch eMsg
               disp(eMsg)
           end
       end     
    end
    end  
end

%update protocol performance measurements
n.packetsToBS(length(n.packetsToBS)+1)=n.packetsToBS(length(n.packetsToBS))-(n.numNodes-n.numDeadNodes);
n.activeNodes(length(n.activeNodes)+1)=n.numNodes-n.numDeadNodes;
sumE=0;
for i=1:n.numNodes
    if (n.node(i).energy)>0
    sumE=sumE+n.node(i).energy;
    end
end
n.consumedEnergy(length(n.consumedEnergy)+1)=(n.initialEnergy*n.numNodes-sumE);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Auxiliary Functions%%%%%%%%%%%%%%%%%%%%%%%
function energy=transmitEnergy(node1,node2,packetlength)
%energy cost (J) to transmit 500 byte message from node1 to node 2
%input example e=transmitEnergyCost(n.node(1),n.BS)

Eelec=50*10^-9; %50 nJ/bit
efs=10*10^-12; %10 pJ/bit/m^2
emp=.0013*10^-12; %pJ/bit/m^4
d0=sqrt(efs/emp);

if nargin<2
    %fixed energy transmission to send Hello to neighbors
    global MAXDIST; % 25 m
    energy=(8*25)*Eelec+(8*25)*efs*MAXDIST^2;
else
    if nargin==3
        len=packetlength;
    else
        len=500;
    end
    %d0 - free space model threshold, 40 m
    %d<do uses freespace , d>d0 multipath
    d=sqrt((node1.x-node2.x)^2+(node1.y-node2.y)^2);
    if d<d0
        energy=(8*len)*Eelec+(8*len)*efs*d^2;
    else
        energy=(8*len)*Eelec+(8*len)*emp*d^4;
    end
end


function energy=receiveEnergy()
%energy cost to transmit from node1 to node 2
%input example e=transmitEnergyCost(n.node(1),n.BS)
%update later

Eelec=50*10^-9; %50 nJ/bit
energy=(8*500)*Eelec;
%energy=50*0.000000001*sqrt((node1.x-node2.x)^2+(node1.y-node2.y)^2);

function d=distance(n1,n2)
d=sqrt((n1.x-n2.x)^2+(n1.y-n2.y)^2);


function energy=aggEnergy(nsignal)
%energy cost (J) to transmit 500 byte message from node1 to node 2
%Data aggrated message at CH
%input example e=transmitEnergyCost(n.node(1),n.BS)

Eda=5*10^-9; %nJ/bit/signal
energy=(8*500)*Eda*nsignal;


function energy=transmitEnergyLEACH(node1,node2,messageType)
%bandwidth=1; %Mb/s
dataMessage=500; %bytes
packetHeader=25; %bytes
if nargin==0
    messageType=1;
end
if messageType==0 %data
    l=(dataMessage+packetHeader)*8; %bits
elseif messageType==1 %logistical stuff
    l=(packetHeader)*8;
end
E_elec=50*10^-9; % J/bit
epsilon_fs=10*10^-12; %J/bit/m^2
epsilon_mp=0.0013*10^-12; %J/bit/m^4
d_o=sqrt(epsilon_fs/epsilon_mp);
if nargin==0
    d=200*2^.5; %guess?? how much energy for broadcast? max dist?
else
    d=sqrt((node1.x-node2.x)^2+(node1.y-node2.y)^2);
end
if d<d_o
	% free space (fs) model
	E_Tx=l*E_elec+l*epsilon_fs*d^2;
else
	% multipath (mp) model
	E_Tx=l*E_elec+l*epsilon_mp*d^4;
end
energy=E_Tx;

function energy=receiveEnergyLEACH(messageType)
dataMessage=500; %bytes
packetHeader=25; %bytes
if messageType==0 %data
    l=(dataMessage+packetHeader)*8; %bits
elseif messageType==1 %for TDMA
    l=(packetHeader)*8;
end
E_elec=50*10^-9; % J/bit
E_Rx=l*E_elec;
energy=E_Rx;

function energy=aggregateDataEnergyLEACH(numNeighbors)
dataMessage=500; %bytes
packetHeader=25; %bytes
l=(dataMessage+packetHeader)*8; %bits
E_DA=5*10^-9; % J/bit/signal
energy=l*E_DA*(numNeighbors+1); %aggregating data from all members of the cluster

