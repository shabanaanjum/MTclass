
%==============
% function oneRoundWSN returns:
% i) Layout of the generated WSN. 
% ii) Graph with the number of messages processed during one round of
% the data gathering. Data are unicasted by each node to the base station(receiver).
% grid=1 (grid topology), grid=0 (random topology)
% receiver: ID of base station
%==============
function oneRoundWSN(numNodes,grid,receiver)
fieldX=500;
fieldY=300;
%Parameters for grid topology
numNodesXY=round(sqrt(numNodes));
step=10;

%% =============Main================
R=calc_R(numNodes,fieldX,fieldY)
[netM,RxTxM]=create_netM(numNodesXY,step,grid,fieldX,fieldY);
figure('Color','w','Position',[100 100 800 500]);
E=printNet(R,netM,fieldX,fieldY);
pause;
RxTxM=simulation(netM,R,fieldX,fieldY,receiver,RxTxM,E,numNodes)
print_RxTxM(RxTxM,numNodes)


%% =============Functions================
function R=calc_R(numNodes,fieldX,fieldY)
        p=sqrt((fieldX^2)+(fieldY^2));
        R=p*sqrt(log10(numNodes)/numNodes);

function [netM,RxTxM]=create_netM(numNodesXY,step,grid,fieldX,fieldY)
ID=1;
for i=1:numNodesXY
    for j=1:numNodesXY
        netM(1,ID)=ID;% inicializaec topologie
        RxTxM(1,ID)=ID; % inicializace matice RxTxM
        if grid==1
            x=step*j+50;
            y=step*i+50;
        else
            x=rand*fieldX;
            y=rand*fieldY;
        end
        netM(2,ID)=x;
        netM(3,ID)=y;
        RxTxM(2,ID)=0;
        RxTxM(3,ID)=0;
        ID=ID+1;
        
    end
end
% Neigbor matrix creation
function E=printNet(R,netM,fieldX,fieldY)
    set(gca,'FontSize',8,'YGrid','off')
    xlabel('\it x \rm [m] \rightarrow')
    ylabel('\it y \rm [m] \rightarrow')
    plot(netM(2,:),netM(3,:),'ko','MarkerSize',5,'MarkerFaceColor','k');
    axis([0 fieldX 0 fieldY]);
    hold all;
    radek=1;
    for j=1:numel(netM(1,:))
        for jTemp=1:numel(netM(1,:))
         X1=netM(2,j);
         Y1=netM(3,j);
         X2=netM(2,jTemp);
         Y2=netM(3,jTemp);
         xSide=abs(X2-X1);
         ySide=abs(Y2-Y1);
         d=sqrt(xSide^2+ySide^2);
         if (d<R)&&(j~=jTemp)
             vertice1=[X1,X2];
             vertice2=[Y1,Y2];
             plot(vertice1,vertice2,'-.k','LineWidth',0.1);
             hold all;
             E(radek,1)=j;
             E(radek,2)=jTemp;
             E(radek,3)=d;
             radek=radek+1;
         end
        end
    end
    v=netM(1,:);
    vv=v';
    s=int2str(vv);
    text(netM(2,:)+1,netM(3,:)+3,s,'FontSize',8,'VerticalAlignment','Baseline');
    hold all;
%% Shortest Path Establishment
function sp=shortestPath(E,sender,receiver)             
        
        [~,sp]=grShortPath(E,sender,receiver);

%% Node highliting
function mark_Ref_Nodes (net,sp,barva)
         for j=1:numel(sp)
             node=sp(j);
             n_X=net(2,node);
             n_Y=net(3,node);
             plot (n_X,n_Y,'bo','LineWidth',3 ,'MarkerEdgeColor', barva,'MarkerSize',6);
         end
            
%% Records unicast transmission into the TxRx matrix
function [RxTxM,sp]=unicast(E,RxTxM,sender,receiver)
        sp=shortestPath(E,sender,receiver);
        for j=1:numel(sp)
            node=sp(j);
            if j==1
                RxTxM(3,node)=RxTxM(3,node)+1;
            elseif j==numel(sp)
                RxTxM(2,node)=RxTxM(2,node)+1;
            else
            RxTxM(2,node)=RxTxM(2,node)+1;
            RxTxM(3,node)=RxTxM(3,node)+1;
            end
               
        end
%% Prints number of processed messages per each node
function print_RxTxM(RxTxM,numNodes)
        figure('Color','w','Position',[100 100 800 500])
        RxTxM(4,:)=RxTxM(2,:)+RxTxM(3,:);
        bar(RxTxM(4,:));
        set(gca,'FontSize',6,'YGrid','off','YGrid','on','XLim',[0 numNodes],'XMinorTick','on');
        xlabel('\it node ID \rm [-] \rightarrow');
        ylabel('\it Number of messages \rm [-] \rightarrow');
        
        hold on;

%% Simulates one round of data gathering    
function RxTxM=simulation(net,R,fieldX,fieldY,receiver,RxTxM,E,numNodes)
        for j=1:numNodes
            sender=j;
            if sender~=receiver
                [RxTxM,sp]=unicast(E,RxTxM,sender,receiver);
                printNet(R,net,fieldX,fieldY);
                mark_Ref_Nodes (net,sp,'r');
                pause(0.01);
                hold off;
            end
        end
    
    
    
    
    
    
    
    
    
    
    
    


