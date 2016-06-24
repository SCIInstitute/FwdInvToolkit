function [x_WBPR_forward,x_WBPR_backward]=WBPR_main(A,y_n,fac,pts,neighb_array,time_interval,pacing_node,last_act_node,QRS_len,ref)
% [x_WBPR_forward,x_WBPR_backward,x_tikh]=WBPR(A,y_n,time_interval,pacing_node,last_act_node,QRS_len)
% This function applies the WBPR method
% A: Forward matrix
% y_n : torso measurements
% time_interval : an array (like m:n) that specifies the QRS time instants  
% pacing_node: the node number of the mesh that is paced
% last_act_node : the last node that is being activated (approximately)
% QRS_len : The number of time instants in the QRS interval (used to
% calculate the reference drift)
% ref: This value specifies the reference drift. In this implementation the
% reference drift is given as an array and the program does not calculate
% it from the solution. The code for calculation of the drift is commented.
% pts: points geometry array N by 3 (without heather)
% fac: face array  M by 3 (without heather)
% neighb_array: nbours file array (without heather)

%
%  For more information, please see: http://software.sci.utah.edu
% 
%  The MIT License
% 
%  Copyright (c) 2009 Scientific Computing and Imaging Institute,
%  University of Utah.
% 
%  
%  Permission is hereby granted, free of charge, to any person obtaining a
%  copy of this software and associated documentation files (the "Software"),
%  to deal in the Software without restriction, including without limitation
%  the rights to use, copy, modify, merge, publish, distribute, sublicense,
%  and/or sell copies of the Software, and to permit persons to whom the
%  Software is furnished to do so, subject to the following conditions:
% 
%  The above copyright notice and this permission notice shall be included
%  in all copies or substantial portions of the Software.
% 
%  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
%  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
%  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
%  DEALINGS IN THE SOFTWARE.
%


% utah epi
%top_nodes=[89:103 283:288 332:336 381:384 428:432 475:480 483 486 488:490];
% T_nij epi
%top_nodes=[1:25 223:231 235:238];
top_nodes = [];
% top_nodes specifies the top nodes of the mesh that we want to stop the
% activation wavefront from propagation

LAMBDA=0.05;
POT_FUNC_DECAY_FACTOR=.3;
WAVE_FRONT_HEIGHT=30;
%
% FIX!!! for different heart geometries
%
neighb_set=find_n_layer_neighb(pacing_node, neighb_array, 2);
initial_pot=zeros(size(pts,1),1);
initial_pot([neighb_set])=-WAVE_FRONT_HEIGHT; % Negative potential of the activated region

[curve_points,curve_start_ptr,No_Contour_Flag]=get_wave_front_new(fac,pts,initial_pot+WAVE_FRONT_HEIGHT/2);

%[interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_points);
  
F=initial_pot;
x_reg=[];
% ref(1)=0;
predicted_pot=[];
 

[u,s,v]=csvd(A);
 
% Forward WBPR
lam = 0;
lam_vector = zeros(length(time_interval),1);

in_curve_flag=2*double(initial_pot>-WAVE_FRONT_HEIGHT/2)-1;

for m=1:length(time_interval)

%    dist=find_dist_from_wave_front_points(interp_sampled_curve_points,pts);

dist=find_dist_from_wave_front_points(curve_points,pts);
    % in_curve_flag=find_in_out_curve(interp_sampled_curve_points,pts);


    % Potential function based on the distance to the wavefront. Specific
    % for the case that negative potential of the activated region is -20

    %  Target geometry size
    target_pts = size(A,2);

    F=WAVE_FRONT_HEIGHT/2*(ones(1,target_pts)-exp(-POT_FUNC_DECAY_FACTOR*dist))'.*in_curve_flag-WAVE_FRONT_HEIGHT/2+ref(m);

    F=[F;zeros(size(A,2)-target_pts,1)]; % zeros added for the case we have a top for the heart mesh
    predicted_pot(:,m)=F;

    tmp=lam;
    %lam=l_curve(u,s,y_n(:,m)-A*F);
    %if(m==57)
    %    lam = tmp;
    %end
    lam_vector(m)=lam;
    lam=LAMBDA;   % fixed lambda was used here instead of l-curve
    x_reg(:,m)=(A'*A+lam^2*eye(size(A,2)))\(A'*y_n(:,time_interval(m))+lam^2*F);
    

    %      if m<2*QRS_len/3
    %          ref(m+1)=x_reg(last_act_node,m);
    %      else
    %          ref(m+1)=20+mean(x_reg(neighb_set,m));
    %      end

    new_thr=-.4*WAVE_FRONT_HEIGHT+ref(m+1);

    activated_top_nodes=find(x_reg(top_nodes,m)>new_thr);
    x_reg(top_nodes(activated_top_nodes),m)=new_thr+.1;


    [curve_points,curve_start_ptr,No_Contour_Flag]=get_wave_front_new(fac,pts,x_reg(:,m)-new_thr);

    in_curve_flag=2*double(x_reg(:,m)>new_thr)-1;


%     if(No_Contour_Flag == 1)
%         x_reg(:,m)=x_reg(:,m-1);
%         continue;
%     end
% 
%     curve_num=length(curve_start_ptr);
%     curve_start_ptr(curve_num+1)=size(curve_points,1)+1;
%     clear curve_cell;
%     for curve_cnt=1:curve_num
%         curve_cell(curve_cnt)={curve_points(curve_start_ptr(curve_cnt):curve_start_ptr(curve_cnt+1)-1,:)};
%         curve_len(curve_cnt)=length(curve_cell{curve_cnt});
%     end
% 
%     [sort_len,ind]=sort(curve_len,'descend');
% 
%     if m<16
%         [interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(1)});
%         if curve_num>1
%             [temp,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(2)});
%             interp_sampled_curve_points=[interp_sampled_curve_points;temp];
%         end
%     else
%         [interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(1)});
%     end
%     [m new_thr curve_start_ptr];
end
x_WBPR_forward=x_reg;

 % Backward WBPR
 
 x_reg=[];
 for m=length(time_interval):-1:1
     
     %dist=find_dist_from_wave_front_points(interp_sampled_curve_points,pts);
     dist=find_dist_from_wave_front_points(curve_points,pts);
     %in_curve_flag=find_in_out_curve(interp_sampled_curve_points,pts);
     
     
     F=WAVE_FRONT_HEIGHT/2*(ones(1,target_pts)-exp(-POT_FUNC_DECAY_FACTOR*dist))'.*in_curve_flag-WAVE_FRONT_HEIGHT/2+ref(m); 
%     F=[F;zeros(size(A,2)-target_pts,1)];
 %    predicted_pot(:,m)=F;
     
     %lam=l_curve(u,s,y_n(:,m)-A*F);
      
     lam=LAMBDA;
     
     x_reg(:,m)=(A'*A+lam^2*eye(size(A,2)))\(A'*y_n(:,time_interval(m))+lam^2*F);

     
%      if m<2*QRS_len/3
%          ref(m+1)=x_reg(last_act_node,m);
%      else
%          ref(m+1)=20+mean(x_reg(neighb_set,m));
%      end

     new_thr=-.6*WAVE_FRONT_HEIGHT+ref(m+1);
     
     activated_top_nodes=find(x_reg(top_nodes,m)<new_thr);
     x_reg(top_nodes(activated_top_nodes),m)=new_thr+.1;
     
     [curve_points,curve_start_ptr,No_Contour_Flag]=get_wave_front_new(fac,pts,x_reg(:,m)-new_thr);

     in_curve_flag=2*double(x_reg(:,m)>new_thr)-1;
     
%      if(No_Contour_Flag == 1)
%          x_WBPR_backward(:,m)=x_WBPR_forward(:,m);
%          m
%          continue;
%      end
%      
%      curve_num=length(curve_start_ptr);
%      curve_start_ptr(curve_num+1)=size(curve_points,1)+1;
%      clear curve_cell;
%      for curve_cnt=1:curve_num
%          curve_cell(curve_cnt)={curve_points(curve_start_ptr(curve_cnt):curve_start_ptr(curve_cnt+1)-1,:)};
%          curve_len(curve_cnt)=length(curve_cell{curve_cnt});
%      end
%          
%      [sort_len,ind]=sort(curve_len,'descend');
%      
%       if m<16
%          [interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(1)});
%          if curve_num>1
%              [temp,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(2)});
%             interp_sampled_curve_points=[interp_sampled_curve_points;temp];
%          end
%       else
%          [interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(1)});
%       end
% 
%      [interp_sampled_curve_points,sampled_curve_points,curve_len]=resample_curve(curve_cell{ind(1)});

     %[m new_thr curve_start_ptr]      
 end
    
x_WBPR_backward=x_reg;

 
