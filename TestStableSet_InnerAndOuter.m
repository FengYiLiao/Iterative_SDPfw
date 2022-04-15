%%TestStableSet both direction

clear;clc;
filename = 'n30p2'
Adj = readmatrix(['SedumiData\StableSet\',filename,'.txt']);%adjacency matrix 
Inner = [];  Outer = [];
Marker = ['o','d','p'];
Color = ['r','g','b'];
line_clr = [0.000000000000000   0.447058826684952   0.741176486015320;
            0.850980401039124   0.325490206480026   0.098039217293262;
            0.000000000000000   0.498039215803146   0.000000000000000;
            0.494117647409439   0.184313729405403   0.556862771511078];
for dx = [1,2,5]
    [In,Out,sdp] = StableSet(Adj,dx);
    Inner = [Inner;In];
    Outer = [Outer;Out];
end

% prob1 = convert_sedumi2mosek(At_sdp,b_sdp,c_sdp,K_sdp);
% [rcode1, res1] = mosekopt('minimize info', prob1);
% OBJ_sdp =  -res1.sol.itr.pobjval;

for num = 1 :3
   if num == 1
       leg1 = plot(1:11,Inner(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
       hold on
       plot(1:11,Outer(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
   elseif num == 2
       leg2 = plot(1:11,Inner(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
       hold on
       plot(1:11,Outer(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
   else
       leg3 = plot(1:11,Inner(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
       hold on
       plot(1:11,Outer(num,:),'Color',line_clr(num,:),'Marker',Marker(num));
   end
end
leg4 = plot(1:11,sdp*ones(11,1),'Color',line_clr(4,:));


%set(gca,'FontSize',20)
xlim([1,11]);
set(gca, 'XTick',  1:1:11)

xlabel('Iteration','Interpreter','latex','FontSize',18);
ylabel('Cost','Interpreter','latex','FontSize',18);
set(gca,'TickLabelInterpreter','latex','fontsize',18);
set(gcf,'Position',[0.18,0.18,0.7,0.7]);
leg = legend([leg4,leg1,leg2,leg3],'$PSD$','$SDD$','$\mathcal{FW}^{10}_{\alpha,2}$',...
    '$\mathcal{FW}^{10}_{\beta,2}$',...
    'Interpreter','latex','Location','northeast','FontSize',20);
set(gcf, 'OuterPosition',  [400 400 600 600])