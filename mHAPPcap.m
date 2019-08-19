function [X,T,P,Qp,Qm,Y,fval,exitflag,status]=mHAPPcap(Linkx,Linkt,Activities,TransNodes,Coeff,LagrM,AffLinks)
%mHAPP is multimodal HAPP (6-24-14) – this script sets up the MILP and solves it using CPLEX
%not implementing activity selection here
%mHAPPcap accounts for Lagrange multiplier for a capacity at a link
%model details can be found at: Chow, J.Y.J., Djavadian, S., 2015. Activity-based market equilibrium for capacitated multimodal transport systems, Transportation Research Part C, Special Issue from ISTTT 21, in press, doi:10.1016/j.trc.2015.04.028.
%Joseph Chow, Ryerson University, 8-12-15
%Inputs - ignoring travel costs
%Linkx -- due to computational cost of complete graph approach, we
%create a link incidence matrix for X with each row being a link ID, and
%cols: r, s, u, w, v at u, v at w, c, t
%Linkt -- similarly, this list is for indexing T,Q variables: cols: r,s,u,v
%Activities is the set of activities to be performed, along with their
%other information: (Pplusx 4 columns)
% Column 1: duration of activity, (d_u)
% Column 2: desired arrival time, (g_u)
% Column 3: early hard time window, (a_u)
% Column 4: late hard time window, (b_u)
%TransNodes is a matrix of multimodal facility information needed to get from one activity to another -- if they are not used, the
%travel is assumed to be walking
% Col 1: mode (v)
% If parking facility, (v=1)
% Col 2: parking fee per minute (p1)
% Col 3: fixed parking fee (p2)
% If fixed schedule transit facility (v=2+)
% Col 2: headway (for uniform dist., 2*d_urs)
% Col 3: scheduled time (a_u=b_u)
%Coeff is a (5+2*Pplus+V)x1 vector of the weights of the objectives,
%for coefficient of penalties, go by early: node 1, 2, ..., then late:
%node 1,2,..
%I will only implement frequency based
%LagrM -- input Lagrange multipliers and associated time: 1x3, where the first col is
%multiplier, and second col is start time, third col is end time
%AffLinks -- index in numX of links affected; e.g. [1; 4; 15] means
%links 1, 4, and 15 are affected
%Outputs
%X - route, N*N
%T - arrival time, N (home start goes first)
%P - deviation from desired arrival, Pplus*2
%Qp - time that a vehicle is picked up at node u from activity r to
%activity s
%Qm - time that a vehicle is dropped off
%Ya, Yb, Yc -- used to capture Lagrange multiplier cost
bigM=10000;
N=size(Activities,1)*2+2;
L=size(TransNodes,1);
numV=max(TransNodes(:,1));
numX=size(Linkx,1);
numtnodes=size(Linkt,1);
Nplus=(N-2)/2; %note that Pminus = Pplus (since this is size of each set)
%objective is all the X's (activity to activity, activity origin, activity
%destination, trans to trans), and then T's (in activities, then trans)
%then T_0 by vehicle and T_2N+1 by vehicle: X,T,P,Qp,Qm
numY=size(AffLinks,1);
numVar=numX+N+numtnodes+Nplus*2+numtnodes*2+numY*3;
lb=zeros(1,numVar);
ub=zeros(1,numVar)+inf;
%upper bound is assumed to be 1440 minutes
ub(1,numX+1:numX+N+numtnodes+Nplus*2+numtnodes*2)=1440;
ub(1,numX+N+numtnodes+Nplus*2+1:numX+N+numtnodes+Nplus*2+numtnodes*2)=0;
for t=1:numtnodes
    for n=1:numX
        if and(Linkx(n,1)==Linkt(t,1),and(Linkx(n,2)==Linkt(t,2),and(Linkx(n,3)==Linkt(t,3),and(Linkx(n,6)==0,and(Linkx(n,4)>0,or(Linkx(n,4)<=Nplus,Linkx(n,4)>N-1))))))
            ub(1,numX+N+numtnodes+Nplus*2+numtnodes+t)=1440;
        end
        if and(Linkx(n,1)==Linkt(t,1),and(Linkx(n,2)==Linkt(t,2),and(Linkx(n,4)==Linkt(t,3),and(Linkx(n,5)==0,or(and(Linkx(n,3)<=Nplus,Linkx(n,3)>0),Linkx(n,3)>N-1)))))
            ub(1,numX+N+numtnodes+Nplus*2+t)=1440;
        end
    end
end
%
%Minimization objective
Obj=zeros(8+numV,numVar);
%Obj z1 - min cost
Obj(1,1:numX)=Coeff(1,1)*Linkx(:,7)';
%Obj z2 - min travel time
Obj(2,1:numX)=Coeff(2,1)*Linkx(:,8)';
%Obj z3 - min parking cost
for i=1:numtnodes
    if Linkt(i,4)==1
        Obj(3,numX+N+numtnodes+Nplus*2+i)=Coeff(3,1)*TransNodes(Linkt(i,3)-N+1,2);
        Obj(3,numX+N+numtnodes+Nplus*2+numtnodes+i)=-Coeff(3,1)*TransNodes(Linkt(i,3)-N+1,2);
    end
end
for j=1:numX
    if Linkx(j,6)==1
        Obj(3,j)=Coeff(3,1)*TransNodes(Linkx(j,4)-N+1,3);
    end
end
%Obj z4 - mode-specific constant
for i=1:numX
    if and(Linkx(i,1)==Linkx(i,3),Linkx(i,6)>0)
        Obj(3+Linkx(i,6),i)=-Coeff(3+Linkx(i,6),1);
    end
end
%Obj z5 - activity return home delay
Obj(4+numV,numX+1+1:numX+1+Nplus)=-Coeff(4+numV,1);
Obj(4+numV,numX+1+Nplus+1:numX+1+2*Nplus)=Coeff(4+numV,1);
%Obj z6 - min length of day
Obj(5+numV,numX+N)=Coeff(5+numV,1);
Obj(5+numV,numX+1)=-Coeff(5+numV,1);
%Obj z7 - early penalty
for u=1:Nplus
    Obj(6+numV,numX+N+numtnodes+u)=Coeff(5+numV+u,1);
    %Obj z8 - late penalty
    Obj(7+numV,numX+N+numtnodes+Nplus+u)=Coeff(5+numV+Nplus+u,1);
end
%Lagrange multiplier
Obj(8+numV,numX+N+numtnodes+Nplus*2+numtnodes*2+numY*2+1:numVar)=LagrM(1,1);
%initiate IP constraint space
Aeq=zeros(N-1+1+(N-1)*(N-1)*L+N-2+L+Nplus+1+Nplus,numVar);
beq=zeros(N-1+1+(N-1)*(N-1)*L+N-2+L+Nplus+1+Nplus,1);
A=zeros(Nplus+2*L*(N-1)*(N-1)+(N-2)*(N-2)+L*(N-1)*(N-2)*2+(N-1)*L*L+N-2+N-2+(N-2)*L*2+(N-1)*(N-2)*2+(N-1)*(N-1)*L*2+numY*3,numVar);
b=zeros(Nplus+2*L*(N-1)*(N-1)+(N-2)*(N-2)+L*(N-1)*(N-2)*2+(N-1)*L*L+N-2+N-2+(N-2)*L*2+(N-1)*(N-2)*2+(N-1)*(N-1)*L*2+numY*3,1);
ineq=0;
ieq=0;
%eq (3) - one flow per activity node
for r=0:N-2
    ieq=ieq+1;
    for n=1:numX
        if and(Linkx(n,1)==r,Linkx(n,3)==r)
            Aeq(ieq,n)=1;
        end
    end
    beq(ieq,1)=1;
end
%eq (4) - return home
for s=1:N-1
    ieq=ieq+1;
    for n=1:numX
        if and(Linkx(n,2)==s,Linkx(n,4)==s)
            Aeq(ieq,n)=1;
        end
    end
    beq(ieq,1)=1;
end
%eq (5) - in equals out
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            ieq=ieq+1;
            count=0;
            for n=1:numX
                if and(Linkx(n,3)==u,and(Linkx(n,1)==r,Linkx(n,2)==s))
                    Aeq(ieq,n)=1;
                    count=1;
                end
                if and(Linkx(n,4)==u,and(Linkx(n,1)==r,Linkx(n,2)==s))
                    Aeq(ieq,n)=-1;
                end
            end
            if count==0
                Aeq(ieq,:)=0;
                ieq=ieq-1;
            end
        end
    end
end
%eq (6) - parked vehicle flow conservation
for u=N:N+L-1
    if TransNodes(u-N+1,1)==1
        ieq=ieq+1;
        count=0;
        for n=1:numX
            if and(Linkx(n,3)==u,Linkx(n,4)>Nplus)
                if or(Linkx(n,4)<=N-1,Linkx(n,6)==1)
                    Aeq(ieq,n)=1;
                    count=1;
                end
            end
            if and(Linkx(n,4)==u,or(and(Linkx(n,3)>=N,Linkx(n,5)==1),or(Linkx(n,3)==0,and(Linkx(n,3)>Nplus,Linkx(n,3)<=N-2)))) %fix 8/25/14
                Aeq(ieq,n)=-1;
            end
        end
        if count==0
            Aeq(ieq,:)=0;
            ieq=ieq-1;
        end
    end
end
% %eq (7) - infeasible actions - implied by the network
%eq (8) - goal arrival
for u=1:Nplus
    ieq=ieq+1;
    Aeq(ieq,numX+1+u)=1;
    Aeq(ieq,numX+N+numtnodes+u)=1;
    Aeq(ieq,numX+N+numtnodes+Nplus+u)=-1;
    beq(ieq,1)=Activities(u,2);
end
%eq (9) - pickup before dropoff
for u=1:Nplus
    ineq=ineq+1;
    A(ineq,numX+1+u)=1;
    A(ineq,numX+1+Nplus+u)=-1;
    b(ineq,1)=-Activities(u,1);
end
%eq (10a) - parking duration
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            if TransNodes(u-N+1,1)==1
                ineq=ineq+1;
                count=0;
                for t=1:numtnodes
                    if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                        A(ineq,numX+N+t)=1;
                        A(ineq,numX+N+numtnodes+Nplus*2+t)=-1;
                    end
                end
                for n=1:numX
                    if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,3)==u))
                        if and(Linkx(n,4)>Nplus,or(Linkx(n,4)<=N-1,Linkx(n,6)==1))
                            A(ineq,n)=bigM;
                            count=1;
                        end
                    end
                end
                if count==0
                    A(ineq,:)=0;
                    ineq=ineq-1;
                else
                    b(ineq,1)=bigM;
                end
            end
        end
    end
end
%eq (10b) - parking duration
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            if TransNodes(u-N+1,1)==1
                ineq=ineq+1;
                count=0;
                for t=1:numtnodes
                    if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                        A(ineq,numX+N+t)=-1;
                        A(ineq,numX+N+numtnodes+Nplus*2+t)=1;
                    end
                end
                for n=1:numX
                    if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,3)==u))
                        if and(Linkx(n,4)>Nplus,or(Linkx(n,4)<=N-1,Linkx(n,6)==1))
                            A(ineq,n)=bigM;
                            count=1;
                        end
                    end
                end
                if count==0
                    A(ineq,:)=0;
                    ineq=ineq-1;
                else
                    b(ineq,1)=bigM;
                end
            end
        end
    end
end
%eq (10c) (originally eq (10b) - parking duration
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            if TransNodes(u-N+1,1)==1
                ineq=ineq+1;
                count=0;
                for t=1:numtnodes
                    if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                        A(ineq,numX+N+t)=-1;
                        A(ineq,numX+N+numtnodes+Nplus*2+numtnodes+t)=1;
                    end
                end
                for n=1:numX
                    if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,4)==u))
                        if or(and(Linkx(n,3)>N-1,Linkx(n,5)==1),or(Linkx(n,3)==0,and(Linkx(n,3)>Nplus,Linkx(n,3)<=N-2)))
                            A(ineq,n)=bigM;
                            count=1;
                        end
                    end
                end
                if count==0
                    A(ineq,:)=0;
                    ineq=ineq-1;
                else
                    b(ineq,1)=bigM;
                end
            end
        end
    end
end
%eq (10d) (originally eq (10b) - parking duration
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            if TransNodes(u-N+1,1)==1
                ineq=ineq+1;
                count=0;
                for t=1:numtnodes
                    if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                        A(ineq,numX+N+t)=1;
                        A(ineq,numX+N+numtnodes+Nplus*2+numtnodes+t)=-1;
                    end
                end
                for n=1:numX
                    if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,4)==u))
                        if or(and(Linkx(n,3)>N-1,Linkx(n,5)==1),or(Linkx(n,3)==0,and(Linkx(n,3)>Nplus,Linkx(n,3)<=N-2)))
                            A(ineq,n)=bigM;
                            count=1;
                        end
                    end
                end
                if count==0
                    A(ineq,:)=0;
                    ineq=ineq-1;
                else
                    b(ineq,1)=bigM;
                end
            end
        end
    end
end
%eq (11a) - time flow
for r=1:N-2
    for s=1:N-2
        ineq=ineq+1;
        count=0;
        A(ineq,numX+1+r)=1;
        A(ineq,numX+1+s)=-1;
        for n=1:numX
            if and(Linkx(n,1)==r,Linkx(n,2)==s)
                if and(Linkx(n,3)==r,Linkx(n,4)==s)
                    A(ineq,n)=bigM;
                    count=1;
                    if r<=Nplus
                        b(ineq,1)=bigM-Activities(r,1)-Linkx(n,8);
                    else
                        b(ineq,1)=bigM-Linkx(n,8);
                    end
                end
            end
        end
        if count==0
            A(ineq,:)=0;
            b(ineq,1)=0;
            ineq=ineq-1;
        end
    end
end
%eq (11b)
for r=1:N-2
    for s=1:N-1
        for w=N:N+L-1
            ineq=ineq+1;
            count=0;
            A(ineq,numX+1+r)=1;
            for t=1:numtnodes
                if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==w))
                    A(ineq,numX+N+t)=-1;
                end
            end
            for n=1:numX
                if and(Linkx(n,1)==r,Linkx(n,2)==s)
                    if and(Linkx(n,3)==r,Linkx(n,4)==w)
                        A(ineq,n)=bigM;
                        count=1;
                        if r<=Nplus
                            b(ineq,1)=bigM-Activities(r,1)-Linkx(n,8);
                        else
                            b(ineq,1)=bigM-Linkx(n,8);
                        end
                    end
                end
            end
            if count==0
                A(ineq,:)=0;
                b(ineq,1)=0;
                ineq=ineq-1;
            end
        end
    end
end
%eq (11c)
for r=0:N-2
    for s=1:N-2
        for u=N:N+L-1
            ineq=ineq+1;
            count=0;
            A(ineq,numX+1+s)=-1;
            for t=1:numtnodes
                if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                    A(ineq,numX+N+t)=1;
                end
            end
            for n=1:numX
                if and(Linkx(n,1)==r,Linkx(n,2)==s)
                    if and(Linkx(n,3)==u,Linkx(n,4)==s)
                        A(ineq,n)=bigM;
                        count=1;
                        if TransNodes(u-N+1,1)>=2
                            b(ineq,1)=bigM-TransNodes(u-N+1,2)/2-Linkx(n,8);
                        else
                            b(ineq,1)=bigM-Linkx(n,8);
                        end
                    end
                end
            end
            if count==0
                A(ineq,:)=0;
                b(ineq,1)=0;
                ineq=ineq-1;
            end
        end
    end
end
%eq (11d)
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            for w=N:N+L-1
                if w~=u
                    ineq=ineq+1;
                    count=0;
                    for t=1:numtnodes
                        if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                            A(ineq,numX+N+t)=1;
                        end
                        if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==w))
                            A(ineq,numX+N+t)=-1;
                        end
                    end
                    for n=1:numX
                        if and(Linkx(n,1)==r,Linkx(n,2)==s)
                            if and(Linkx(n,3)==u,Linkx(n,4)==w)
                                A(ineq,n)=bigM;
                                count=1;
                                if TransNodes(u-N+1,1)>=2
                                    b(ineq,1)=bigM-TransNodes(u-N+1,2)/2-Linkx(n,8);
                                else
                                    b(ineq,1)=bigM-Linkx(n,8);
                                end
                            end
                        end
                    end
                    if count==0
                        A(ineq,:)=0;
                        b(ineq,1)=0;
                        ineq=ineq-1;
                    end
                end
            end
        end
    end
end
%eq (11e) -- from origin
for s=1:N-2
    ineq=ineq+1;
    count=0;
    A(ineq,numX+1)=1;
    A(ineq,numX+1+s)=-1;
    for n=1:numX
        if and(Linkx(n,1)==0,Linkx(n,2)==s)
            if and(Linkx(n,3)==0,Linkx(n,4)==s)
                A(ineq,n)=bigM;
                count=1;
                b(ineq,1)=bigM-Linkx(n,8);
            end
        end
    end
    if count==0
        A(ineq,:)=0;
        b(ineq,1)=0;
        ineq=ineq-1;
    end
end
%eq (11f) - to dest
for r=1:N-2
    ineq=ineq+1;
    count=0;
    A(ineq,numX+1+r)=1;
    A(ineq,numX+N)=-1;
    for n=1:numX
        if and(Linkx(n,1)==r,Linkx(n,2)==N-1)
            if and(Linkx(n,3)==r,Linkx(n,4)==N-1)
                A(ineq,n)=bigM;
                count=1;
                if r<=Nplus
                    b(ineq,1)=bigM-Activities(r,1)-Linkx(n,8);
                else
                    b(ineq,1)=bigM-Linkx(n,8);
                end
            end
        end
    end
    if count==0
        A(ineq,:)=0;
        b(ineq,1)=0;
        ineq=ineq-1;
    end
end
%eq (11g)
for s=1:N-1
    for w=N:N+L-1
        ineq=ineq+1;
        count=0;
        A(ineq,numX+1)=1;
        for t=1:numtnodes
            if and(Linkt(t,1)==0,and(Linkt(t,2)==s,Linkt(t,3)==w))
                A(ineq,numX+N+t)=-1;
            end
        end
        for n=1:numX
            if and(Linkx(n,1)==0,Linkx(n,2)==s)
                if and(Linkx(n,3)==0,Linkx(n,4)==w)
                    A(ineq,n)=bigM;
                    count=1;
                    b(ineq,1)=bigM-Linkx(n,8);
                end
            end
        end
        if count==0
            A(ineq,:)=0;
            b(ineq,1)=0;
            ineq=ineq-1;
        end
    end
end
%eq (11h)
for r=0:N-2
    for u=N:N+L-1
        ineq=ineq+1;
        count=0;
        A(ineq,numX+N)=-1;
        for t=1:numtnodes
            if and(Linkt(t,1)==r,and(Linkt(t,2)==N-1,Linkt(t,3)==u))
                A(ineq,numX+N+t)=1;
            end
        end
        for n=1:numX
            if and(Linkx(n,1)==r,Linkx(n,2)==N-1)
                if and(Linkx(n,3)==u,Linkx(n,4)==N-1)
                    A(ineq,n)=bigM;
                    count=1;
                    if TransNodes(u-N+1,1)>=2
                        b(ineq,1)=bigM-TransNodes(u-N+1,2)/2-Linkx(n,8);
                    else
                        b(ineq,1)=bigM-Linkx(n,8);
                    end
                end
            end
        end
        if count==0
            A(ineq,:)=0;
            b(ineq,1)=0;
            ineq=ineq-1;
        end
    end
end
%eq (12a) - activity time windows
for r=1:Nplus
    ineq=ineq+1;
    A(ineq,numX+1+r)=-1;
    b(ineq,1)=-Activities(r,3);
end
%eq (12b) - activity time windows
for r=1:Nplus
    ineq=ineq+1;
    A(ineq,numX+1+r)=1;
    b(ineq,1)=Activities(r,4);
end
%%%I'm leaving out the hard time windows for transit for now
%eq (13a) - to keep Qp's at 0 if not visited
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            ineq=ineq+1;
            count=0;
            for t=1:numtnodes
                if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                    A(ineq,numX+N+numtnodes+2*Nplus+numtnodes+t)=1;
                end
            end
            for n=1:numX
                if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,3)==u))
                    if Linkx(n,6)==0 %and(Linkx(n,4)<=N-1,)
                        count=1;
                        A(ineq,n)=-bigM;
                    end
                end
            end
            if count==0
                A(ineq,:)=0;
                ineq=ineq-1;
            end
        end
    end
end
%eq (13b) - to keep T's at 0 if not visited
for r=0:N-2
    for s=1:N-1
        for u=N:N+L-1
            ineq=ineq+1;
            count=0;
            for t=1:numtnodes
                if and(Linkt(t,1)==r,and(Linkt(t,2)==s,Linkt(t,3)==u))
                    A(ineq,numX+N+t)=1;
                end
            end
            for n=1:numX
                if and(Linkx(n,1)==r,and(Linkx(n,2)==s,Linkx(n,3)==u))
                    count=1;
                    A(ineq,n)=-bigM;
                end
            end
            if count==0
                A(ineq,:)=0;
                ineq=ineq-1;
            end
        end
    end
end
%new constraint for Y_uw^a
for n=1:numY
    ineq=ineq+1;
    A(ineq,AffLinks(n,1))=bigM;
    A(ineq,numX+N+numtnodes+Nplus*2+numtnodes*2+n)=-bigM;
    if Linkx(AffLinks(n,1),3)<=N-1
        A(ineq,numX+1+Linkx(AffLinks(n,1),3))=1;
    else
        for t=1:numT
            if and(Linkx(AffLinks(n,1),1)==Linkt(t,1),and(Linkx(AffLinks(n,1),2)==Linkt(t,2),Linkx(AffLinks(n,1),3)==Linkt(t,3)))
                A(ineq,numX+N+t)=1;
            end
        end
    end
    b(ineq,1)=LagrM(1,2)-Linkx(AffLinks(n,1),8)+bigM;
end
%new constraint for Y_uw^b
for n=1:numY
    ineq=ineq+1;
    A(ineq,AffLinks(n,1))=bigM;
    A(ineq,numX+N+numtnodes+Nplus*2+numtnodes*2+numY+n)=-bigM;
    if Linkx(AffLinks(n,1),3)<=N-1
        A(ineq,numX+1+Linkx(AffLinks(n,1),3))=-1;
    else
        for t=1:numT
            if and(Linkx(AffLinks(n,1),1)==Linkt(t,1),and(Linkx(AffLinks(n,1),2)==Linkt(t,2),Linkx(AffLinks(n,1),3)==Linkt(t,3)))
                A(ineq,numX+N+t)=-1;
            end
        end
    end
    b(ineq,1)=-LagrM(1,3)+bigM;
end
%new constraint for Y_uw^c
for n=1:numY
    ineq=ineq+1;
    A(ineq,numX+N+numtnodes+Nplus*2+numtnodes*2+n)=1;
    A(ineq,numX+N+numtnodes+Nplus*2+numtnodes*2+numY+n)=1;
    A(ineq,numX+N+numtnodes+Nplus*2+numtnodes*2+numY*2+n)=-1;
    b(ineq,1)=1;
end
Aeq=Aeq(1:ieq,:);
A=A(1:ineq,:);
beq=beq(1:ieq,1);
b=b(1:ineq,1);
%CPLEX MILP Solver
varArray=[];
for i=1:numX
    varArray=[varArray 'B'];
end
for i=numX+1:numX+N+numtnodes+Nplus*2+numtnodes*2
    varArray=[varArray 'C'];
end
for i=numX+N+numtnodes+Nplus*2+numtnodes*2+1:numVar
    varArray=[varArray 'B'];
end
f=sum(Obj,1)';
[Z,fval,exitflag,status]=cplexmilp(f,A,b,Aeq,beq,[],[],[],lb',ub',varArray);
X=Z(1:numX,1);
T=Z(numX+1:numX+N+numtnodes,1);
P=Z(numX+N+numtnodes+1:numX+N+numtnodes+2*Nplus,1);
Qp=Z(numX+N+numtnodes+Nplus*2+1:numX+N+numtnodes+Nplus*2+numtnodes);
Qm=Z(numX+N+numtnodes+Nplus*2+numtnodes+1:numX+N+numtnodes+Nplus*2+2*numtnodes);
Y=Z(numX+N+numtnodes+Nplus*2+2*numtnodes+1:numVar);