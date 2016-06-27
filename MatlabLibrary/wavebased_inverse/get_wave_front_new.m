function [curve_points,curve_start_ptr,No_Contour_Flag]=get_wave_front_new(fac,pts,potval)

%[curve_points,curve_start_ptr,No_Contour_Flag]=get_wave_front_new(fac,pts,potval)
% Returns the zero level contours in curve_points (ordered properly) as an array of points (N by 2)
% curve_start_ptr specifies the start points of the curves  in curve_points
% No_Contour_Flag will be 1 if no contour found


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

potval(find(potval==0))=0.01; %this is to prevent having a cotour on a node
No_Contour_Flag=0;
cnt=1;
for i=1:size(fac,1)
    if (max(potval(fac(i,:)))<0 | min(potval(fac(i,:)))>0)
        continue;
    else
        [max_val,max_ind]=max(potval(fac(i,:)));
        [min_val,min_ind]=min(potval(fac(i,:)));
        [temp,third_ind]=min(ismember(fac(i,:),fac(i,[max_ind min_ind])));
        
        if (potval(fac(i,third_ind))<=0)
            act_pts(cnt,1:3)=(potval(fac(i,max_ind))*pts(fac(i,third_ind),:)-potval(fac(i,third_ind))*pts(fac(i,max_ind),:))./(potval(fac(i,max_ind))-potval(fac(i,third_ind)));
            wave_front_nodes(cnt,1:2)=sort(fac(i,[max_ind third_ind]));
        else
            act_pts(cnt,1:3)=(potval(fac(i,third_ind))*pts(fac(i,min_ind),:)-potval(fac(i,min_ind))*pts(fac(i,third_ind),:))./(potval(fac(i,third_ind))-potval(fac(i,min_ind)));
            wave_front_nodes(cnt,1:2)=sort(fac(i,[min_ind third_ind]));
        end
        act_pts(cnt,4:6)=(potval(fac(i,max_ind))*pts(fac(i,min_ind),:)-potval(fac(i,min_ind))*pts(fac(i,max_ind),:))./(potval(fac(i,max_ind))-potval(fac(i,min_ind)));
        wave_front_nodes(cnt,3:4)=sort(fac(i,[max_ind min_ind]));
        
        cnt=cnt+1;
    end
end


if cnt==1
    No_Contour_Flag = 1;
end

if (No_Contour_Flag == 1)
    curve_points = [];
    curve_start_ptr = [];
    return;
end

curve_points(1,:)=act_pts(1,1:3);
curve_points(2,:)=act_pts(1,4:6);
cnt2=2;
cur_pos=1;
curve_num=1;
curve_start_ptr(1)=1;
cur_nodes=wave_front_nodes(1,3:4);
wave_front_nodes(1,3:4)=[0 0];
over_flow_cnt=0;
while cnt2<cnt-1 && over_flow_cnt<1000
    if (cur_nodes == [0 0])
        curve_num=curve_num+1;
        temp=find(sum(wave_front_nodes')'~=0);
        curve_start_ptr(curve_num)=cnt2-1;
        curve_points(cnt2-1,:)=act_pts(temp(1),1:3);
        curve_points(cnt2,:)=act_pts(temp(1),4:6);
        cur_nodes=wave_front_nodes(temp(1),3:4);
        wave_front_nodes(temp(1),3:4)=[0 0];
    end
        
    [is_mem,loc]=ismember(cur_nodes,wave_front_nodes(:,1:2),'rows');
    if is_mem==0
        [is_mem,loc]=ismember(cur_nodes,wave_front_nodes(:,3:4),'rows');
        if is_mem==1
            cnt2=cnt2+1;
            curve_points(cnt2,:)=act_pts(loc,1:3);
            cur_nodes=wave_front_nodes(loc,1:2);            
            wave_front_nodes(loc,:)=[0 0 0 0];
            
        end
    else
        cnt2=cnt2+1;
        curve_points(cnt2,:)=act_pts(loc,4:6);
        cur_nodes=wave_front_nodes(loc,3:4);
        wave_front_nodes(loc,:)=[0 0 0 0];
    end   
    over_flow_cnt=over_flow_cnt+1;
end
