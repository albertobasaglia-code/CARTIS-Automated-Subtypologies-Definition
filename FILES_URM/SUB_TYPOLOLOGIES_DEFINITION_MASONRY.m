%% AUTOMATED SUBTYPOLOGIES DEFINITION BASED ON CARTIS DATA
% Alberto Basaglia, University "G. D'Annunzio" of Chieti-Pescara (2021)

% UNREINFORCED MASONRY

clear all
close all
clc

%% SUPPORTING FILES

Municipalities=readtable('Municipalities.txt','PreserveVariableNames',true);
% List of Italian Municipalities (as of JAN 2021)

%% INPUT DATA

d=dir('data*.csv');
CARTIS=readtable(char((d.name)),'PreserveVariableNames',true);

% Listing of Municipalities provided (check)
Municipalities_CARTIS=unique(CARTIS(:,7));
% Assignation of ID and ID Province
for j=1:size(Municipalities_CARTIS,1)
    Def(j,1:2)=Municipalities(strcmp(Municipalities.Comune, Municipalities_CARTIS.comune(j)), 2:3);
end
Municipalities_CARTIS=[Municipalities_CARTIS Def];

%% CREATION OF MATRIX FOR DISAGGREGATION

A(:,1)=table2array(CARTIS(:,3));

% Municipalities conversion in numbers (table to double)
for j=1:size(CARTIS,1)
    municip_provinces(j,1:2)=Municipalities(strcmp(Municipalities.Comune, CARTIS.comune(j)), 2:3);
end
A=[A table2array(municip_provinces)];

% Conversion URM typologies in numbers (table to double)
tipology=CARTIS(:,8);
tipology(:,1)=regexprep(tipology{:,1},'','');
tipology=(strcmp(tipology{:, 1},'A1.1 - Ciottoli con tessitura disordinata nel parametro')*1 ...
       +strcmp(tipology{:,1},'A1.2 - Ciottoli con tessitura ordinata nel parametro')*2 ...
       +strcmp(tipology{:,1},'A1.3 - Ciottoli e mattoni')*3 ...
       +strcmp(tipology{:,1},'A1.4 - Ciottoli e mattoni con ricorsi in laterizio')*4 ...
       +strcmp(tipology{:,1},'A2.1 - Pietrame con tessitura disordinata nel parametro')*5 ...
       +strcmp(tipology{:,1},'A2.2 - Pietrame con tessitura ordinata nel parametro')*6 ...
       +strcmp(tipology{:,1},'A2.3 - Muratura disordinata con embrici e calcare')*7 ...
       +strcmp(tipology{:,1},'A2.4 - Pietrame con ricorsi in laterizio')*8+strcmp(tipology{:,1},'B1.1 - Senza Ricorsi')*9 ...
       +strcmp(tipology{:,1},'B1.2 - Con Ricorsi')*10+strcmp(tipology{:,1},'B2.1 - Senza Ricorsi')*11 ...
       +strcmp(tipology{:,1},'B2.2 - Con Ricorsi')*12+strcmp(tipology{:,1},'C1.1 - Senza Ricorsi')*13 ...
       +strcmp(tipology{:,1},'C1.2 - Con Ricorsi')*14+strcmp(tipology{:,1},'C2.0 - Mattoni')*15 ...
       +strcmp(tipology{:,1},'C2.1 - Mattoni in cls')*16);
% NOTE: typo in "parametro" which should be "paramento" (from pgAdmin database)
A=[A tipology];

disp('Parameters for disaggregation (pgAdmin JAN 2021).')
disp(' PAR 1: Number of stories')
disp(' PAR 2: Construction period')
disp(' PAR 3: Ring beams and tie rods')
disp(' PAR 4: Slab type')
disp(' PAR 5: Presence of mixed structures')
disp(' PAR 6: Regularity in plan')
disp(' PAR 7: Regularity in elevation')
disp(' PAR 8: Roof material')
disp(' PAR 9: Roof spread')
disp(' PAR 10: Structural strengthening')
disp(' PAR 11: Vault type')
disp(' PAR 12: Average storey height')
disp(' PAR 13: Average ground floor height')
disp(' PAR 14: Percentage of openings')
disp(' PAR 15: Percentage of openings at the ground floor')
disp(' PAR 16: Average floor area')
disp(' PAR 17: Average wall thickness at the ground floor')
disp(' PAR 18: Average distance between walls')

PAR=input('\nInsert number of selected parameters in [], space delimited:\n');
clc

% Columns   [ 1 ] [ 2 ] [ 3 ] [ 4 ] [ 5 ] [ 6 ] [ 7 ] [ 8 ] [ 9 ] [ 10 ] [ 11 ] [ 12 ] [ 13 ] [ 14 ] [  15  ] [  16  ] [  17  ] [  18  ] 
PAR_columns=[9 20 22 35 37 38 39 47 49 54 56 58 60 62 64 67 69 70 72  75 76  83 86  89 91  94 96 100 102  106 108  123 125  127 129 132];
% NOTE: Second column of PAR 3 is NOT N.d. (data not available) but "Absent"
%       Last column of PAR 10 is considered even if it's N.d. as it will become "None" (see below)
PAR_selected=[];

for j=1:18
    % 12 Parameters (max)
    if nnz(PAR(:)==j)>0
       PAR_selected=[PAR_selected size(A,2)+1];
       % Position of first column of selected parameter in matrix A
       A=[A table2array(CARTIS(:,PAR_columns(2*j-1):PAR_columns(2*j)))];
       PAR_selected=[PAR_selected size(A,2)];
       % Position of last column of selected parameter in matrix A
    end
end

% FIX PARAMETERS
% Fix Parameter 3: Ring beams and tie rods
if nnz(PAR(:)==3)>0
    pos=find(PAR(:)==3);
    pos=2*pos;
    % Position of last column of the selected parameters in A
    for k=1:size(A,1)
        if nnz(A(k,PAR_selected(pos-1):PAR_selected(pos)))==0
        % Both values are equal to zero (i.e. data is missing)
        elseif nnz(A(k,PAR_selected(pos-1):PAR_selected(pos)))>0 & abs(sum(A(k,PAR_selected(pos-1):PAR_selected(pos)))-A(k,1))>0
        % One value greater than zero, fix may be needed
            z=find(A(k,:)==0);
            z=z(1,find(z>=PAR_selected(pos-1) & z<=PAR_selected(pos)));
            % Position of cell with value equal to zero
            A(k,z)=A(k,1)-max(A(k,PAR_selected(pos-1):PAR_selected(pos)));
            % Value is replaced so that sum is equal to the no. of buildings in compartment
        end
    end
end      
% Fix Parameter 10: Structural strengthening
if nnz(PAR(:)==10)>0
    pos=find(PAR(:)==10);
    pos=2*pos;
    % Position of last column of the selected parameters in A
    for k=1:size(A,1)
        if abs(sum(A(k,PAR_selected(pos-1):PAR_selected(pos)))-A(k,1))>0
            % NOTE: N.d. values (last column) have to be considered since the option
            % "None" is not included in pgAdmin
            A(k,PAR_selected(pos))=max(0,A(k,1)-sum(A(k,PAR_selected(pos-1):PAR_selected(pos)-1)));
            % LAST column becomes "No structural strengthening"
            % NOTE: max(0,:) is used to avoid negative values (for splitting no. of buildings in different entries)
        else
            % Data was N.d.
            A(k,PAR_selected(pos))=0;
            % As N.d. is not considered, value is set to zero (e.g. information missing)
        end
    end
end

% MATRIX A
% Column 1: no. Buildings in Compartment
% Column 2: ID Municipality
% Column 3: ID Province
% Column 4: ID Typology
% Columns >4: Selected Parameters

%% MATRIX DISAGGREGATION

B=[];
% Disaggregated Matrix
Error=[];
% Error Matrix (due to data missing or inconsistent)
% Data inconsistency is due to value in parameter(s) different from total no. of Buildings

for j=1:size(Municipalities_CARTIS,1)
    % Cycle for each Municipality in the data provided
    B_temp=A(find(A(:,2)==table2array(Municipalities_CARTIS(j,2))),:);
    % Temporary matrix with all compartments of a single Municipality
    
    for k=1:size(B_temp,1)
        % Cycle for each compartment of the Municipality
        pos=1;
        z=[];
        z_sum=[];
        tol=5;
        % Maximum difference between no. of Buildings in compartment
        % (column 1) and sum of buildings in parameter(s)
        
        while pos<length(PAR_selected)
            z=[z nnz(B_temp(k,PAR_selected(pos):PAR_selected(pos+1)))];
            z_sum=[z_sum sum(B_temp(k,PAR_selected(pos):PAR_selected(pos+1)))];
            pos=pos+2;
        end
        
        if nnz(not(z))>0 && B_temp(k,4)>0
            % CASE 1: Data is missing           
            Error=[Error; [Municipalities_CARTIS(j,1) array2table(k) array2table(B_temp(k,1)) {'Data missing'}]];
        elseif nnz(not(z))>0 && B_temp(k,4)==0
            % CASE 2: Data and typology is missing           
            Error=[Error; [Municipalities_CARTIS(j,1) array2table(k) array2table(B_temp(k,1)) {'Data and typology missing'}]];            
        elseif max(abs(z_sum-B_temp(k,1)))>tol && B_temp(k,4)>0
            % CASE 3: Data inconsistency 
            Error=[Error; [Municipalities_CARTIS(j,1) array2table(k) array2table(B_temp(k,1)) {'Data inconsistent'}]];
        elseif max(abs(z_sum-B_temp(k,1)))>tol && B_temp(k,4)==0
            % CASE 4: Data inconsistency and typology missing
            Error=[Error; [Municipalities_CARTIS(j,1) array2table(k) array2table(B_temp(k,1)) {'Data inconsistent and typology missing'}]];
        elseif B_temp(k,4)==0
            % CASE 5: Data inconsistency and typology missing
            Error=[Error; [Municipalities_CARTIS(j,1) array2table(k) array2table(B_temp(k,1)) {'Typology missing'}]];
        else
            % CASE 6: Data is present and consistent: proceed with disaggregation
            B_dis=zeros(prod(z),size(B_temp,2)-1);
            % Temporary matrix with disaggregated compartment
            % ID Province are dropped
            B_dis(:,2)=(B_temp(k,2));
            % ID Municipality is copied
            B_dis(:,3)=(B_temp(k,4));
            % ID Typology is copied
            znot=find(B_temp(k,:));
            
            % Defining values for disaggregation
            % NOTE: Additional cycles must be added if database is updated
            w=1;
            if nnz(PAR(:)==1)>0
                z1=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
               z1=0;
            end
            if nnz(PAR(:)==2)>0
                z2=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
            z2=0;
            end
            if nnz(PAR(:)==3)>0
                z3=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
               w=w+2;
            else
                z3=0;
            end        
            if nnz(PAR(:)==4)>0
                z4=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z4=0;
            end        
            if nnz(PAR(:)==5)>0
                z5=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z5=0;
            end
            if nnz(PAR(:)==6)>0
                z6=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z6=0;
            end        
            if nnz(PAR(:)==7)>0
                z7=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
               z7=0;
            end        
            if nnz(PAR(:)==8)>0
               z8=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z8=0;
            end        
            if nnz(PAR(:)==9)>0
                z9=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z9=0;
            end        
            if nnz(PAR(:)==10)>0
                z10=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z10=0;
            end        
             if nnz(PAR(:)==11)>0
                z11=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z11=0;
            end       
            if nnz(PAR(:)==12)>0
                z12=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z12=0;
            end           
            if nnz(PAR(:)==13)>0
                z13=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z13=0;
            end            
            if nnz(PAR(:)==14)>0
                z14=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z14=0;
            end            
            if nnz(PAR(:)==15)>0
                z15=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z15=0;
            end            
            if nnz(PAR(:)==16)>0
                z16=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z16=0;
            end            
            if nnz(PAR(:)==17)>0
                z17=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z17=0;
            end            
            if nnz(PAR(:)==18)>0
                z18=znot(znot>=PAR_selected(w) & znot<=PAR_selected(w+1));
                w=w+2;
            else
                z18=0;
            end 
            comb=allcomb(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18);
            comb(:,~any(comb,1))=[];
            % Columns with 0 are removed
            comb1=comb-1;
            % As ID Province is dropped
        
            for l=1:prod(z)
                B_dis(l,comb1(l,:))=1;
                % A "1" is put in designated columns (see comb1) to mark the value
                % of selected parameters
                B_dis(l,1)=ceil(B_temp(k,1)*prod(B_temp(k,comb(l,:))./B_temp(k,1)));
                % No. of Buildings in disaggregated Compartment
                % (equally distributed according to distribution among selected parameter)
                % comb is used instead of comb1 as it refers to B_temp (with ID Municipality and ID Province)
            end
        
            B=[B; B_dis];
            % Output    
            
        end
       
    end
    % Cycle repeats for next Municipality
end
PAR_selected=PAR_selected-1;
% As ID Province is dropped


Error.Properties.VariableNames={'Municipality','Compartment','No. of Buildings','Error Type'};

% DISPLAY WARNING (if any)
if isempty(B)==1
    % Disaggregated matrix is empty
    disp('It was not possible to proceed with the disaggregation. Please check "Error" for more details.\n\n')
    return
    % process is halted
elseif sum(table2array(Error(:,3)))/sum(B(:,1))>0.5
    % Less than 50% of buildings were disaggregated
    disp('Warning: less than 50% of buildings were disaggregated.')
    check=input('\nDo you wish to continue? (Yes/No)\n','s');
    if strcmpi(check,'Yes')
        fprintf('\nThe disaggregation is complete. %d buildings were considered (%.2f%% of the total).\n\n',sum(B(:,1)),min(100,(sum(B(:,1))/sum(table2array(CARTIS(:,3))))*100));
    elseif strcmpi(check,'No')
        fprintf('\nThe process has been halted by the user. Please check "Error" for more details.\n\n')
        return
        % process is halted
    end     
else
    fprintf('The disaggregation is complete. %d buildings were considered (%.2f%% of the total).\n\n',sum(B(:,1)),(sum(B(:,1))/sum(table2array(CARTIS(:,3))))*100);
end

% MATRIX B
% Column 1: no. Buildings in Compartment
% Column 2: ID Municipality
% Column 3: ID Typology
% Columns >3: Selected Parameters

%% PARAMETERS' AGGREGATION (OPTIONAL)

agg=input('Perform parameters aggregation? (Yes/No)\n','s');
if strcmpi(agg,'Yes')
    agg=1;
elseif strcmpi(agg,'No')
    agg=0;
end

% Performing aggregation
if agg==1        
    % Typology
    % A1.1, A1.2, A2.1, A2.2 = Irregular masonry - Without courses
    tip1=find(B(:,3)==1 | B(:,3)==2 | B(:,3)==5 | B(:,3)==6);
    % A1.3, A1.4, A2.3, A2.4 = Irregular masonry - With courses
    tip2=find(B(:,3)==3 | B(:,3)==4 | B(:,3)==7 | B(:,3)==8); 
    % B1.1, B2.1 = Rough-hewn masonry - Without courses
    tip3=find(B(:,3)==9 | B(:,3)==11);  
    % B1.2, B2.2 = Rough-hewn masonry - With courses
    tip4=find(B(:,3)==10 | B(:,3)==12); 
    % C1.1 = Square cut stone - Without courses
    tip5=find(B(:,3)==13);
    % C1.2 = Square cut stone - With courses
    tip6=find(B(:,3)==14);
    % C2.0, C2.1 = Blocks
    tip7=find(B(:,3)==15 | B(:,3)==16);
    % Assignation of new typologies
    B(tip1,3)=1;
    B(tip2,3)=2;
    B(tip3,3)=3;
    B(tip4,3)=4;
    B(tip5,3)=5;
    B(tip6,3)=6;
    B(tip7,3)=7;
    
    % Parameter 2: Construction period
    if nnz(PAR==2)>0
        pos=PAR_selected(2*find(PAR==2)-1);
        % Find position of first column of Parameter 2 in B
        B=[B(:,1:pos+2) sum(B(:,pos+3:pos+6),2) sum(B(:,pos+7:pos+13),2) B(:,pos+14:end)];
        % Construction periods -> After 1945, 2 periods
        %                         1946-1981 and Post-1981 (evolution of codes)
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-9;
        % Updating of PAR_selected given the condensed number of columns
    end
    
    % Parameter 4: Slab type
    if nnz(PAR==4)>0
        pos=PAR_selected(2*find(PAR==4)-1);
        % Find position of first column of Parameter 4 in B
        B=[B(:,1:pos-1) sum(B(:,pos:pos+2),2) sum(B(:,pos+3:pos+5),2) sum(B(:,pos+6:pos+8),2) B(:,pos+9:end)];
        % Slab type -> Flexible, Semi-rigid, Rigid
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-6;
        % Updating of PAR_selected given the condensed number of columns
    end
    
    % Parameter 8: Roof material
    if nnz(PAR==8)>0
        pos=PAR_selected(2*find(PAR==8)-1);
        % Find position of first column of Parameter 8 in B 
        B=[B(:,1:pos-1) sum(B(:,pos:pos+1),2) sum(B(:,pos+2:pos+3),2) B(:,pos+4:end)];
        % Roof material -> Wood, steel = Light
        %                  RC, Masonry = Heavy
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-2;
        % Updating of PAR_selected given the condensed number of columns
    end    
    
    % Parameter 17: Average wall thickness at the ground floor
    if nnz(PAR==17)>0
        pos=PAR_selected(2*find(PAR==17)-1);
        % Find position of first column of Parameter 17 in B 
        B=[B(:,1:pos-1) sum(B(:,pos:pos+1),2) B(:,pos+2:end)];
        % Wall thickness -> < 40 cm, 40-60 cm = < 60 cm
        %                   > 60 cm
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-1;
        % Updating of PAR_selected given the condensed number of columns
    end    
    
    % Parameter 18: Average distance between walls
    if nnz(PAR==18)>0
        pos=PAR_selected(2*find(PAR==18)-1);
        % Find position of first column of Parameter 18 in B
        B=[B(:,1:pos-1) sum(B(:,pos:pos+1),2) sum(B(:,pos+2:pos+3),2) B(:,pos+4:end)];
        % Distance between walls -> < 4 m, 4-5 m = < 5 m
        %                           5-6 m, > 6 m = > 6 m
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-2;
        % Updating of PAR_selected given the condensed number of columns
    end
end

%% SUMMING EQUAL ROWS

C=[];
% Output matrix (double)

for j=1:max(B(:,3))
    % Cycle for each building Typology
    C1=[];
    % Temporary matrix of all buildings with the same typology
    C1=find(B(:,3)==j);
    C1=B(C1,[1,4:end]);
    while isempty(C1)==0
        C1diff=C1-[0 C1(1,2:end)];
        % To look for other rows with buildings having the same
        % characteristics of the first 
        % NOTE: first column is not considered, i.e. no. of buildings
        d1=find(all(C1diff(:,2:end)==0,2));
        % Rows of buildings with the same characteristics
        C=[C; [sum(C1diff(d1,1)) j C1(1,2:end)]];
        % ID Municipality is dropped
        C1(d1,:)=[];
        % Removing rows summed
    end
end
PAR_selected=PAR_selected-1;
% As ID Municipality is dropped

% Sorting and adding percentages
C=[C (C(:,1)./sum(C(:,1)))*100];
C=sortrows(C,size(C,2),'descend');
C=[C cumsum(C(:,end))];

TOT=sum(C(:,1));

% MATRIX C
% Column 1: no. Buildings in Compartment
% Column 2: ID Typology
% Columns >2: Selected Parameters
% Column end-1: Percentage on dataset
% Column end: Cumulative percentage on dataset

%% OUTPUT

Output=[];

% NO. OF BUILDINGS AND TYPOLOGY
Output=[Output array2table(C(:,1))];
Output.Properties.VariableNames={'No. of Builidngs'};

tipology=[];
if agg==0
    % Without aggregation
    for j=1:size(C,1)
        if C(j,2)==1
            tipology=[tipology; cell2table({'A1.1 Pebble wall with irregular texture'})];
        elseif C(j,2)==2
            tipology=[tipology; cell2table({'A1.2 Pebble wall with regular texture'})];
        elseif C(j,2)==3
            tipology=[tipology; cell2table({'A1.3 Brick and pebble wall'})];
        elseif C(j,2)==4
            tipology=[tipology; cell2table({'A1.4 Brick and pebble wall with brick courses'})];
        elseif C(j,2)==5
            tipology=[tipology; cell2table({'A2.1 Stone wall with irregular texture'})];
        elseif C(j,2)==6
            tipology=[tipology; cell2table({'A2.2 Stone wall with regular texture'})];
        elseif C(j,2)==7
            tipology=[tipology; cell2table({'A2.3 Irregular masonry with brick tiles and limestone'})];
        elseif C(j,2)==8
            tipology=[tipology; cell2table({'A2.4 Stone wall with brick courses'})];
        elseif C(j,2)==9
            tipology=[tipology; cell2table({'B1.1 Stacked slab stone wall - without courses'})];
        elseif C(j,2)==10
            tipology=[tipology; cell2table({'B1.2 Stacked slab stone wall - with courses'})];
        elseif C(j,2)==11
            tipology=[tipology; cell2table({'B2.1 Stone wall with pseudo regular blocks - without courses'})];
        elseif C(j,2)==12
            tipology=[tipology; cell2table({'B2.2 Stone wall with pseudo regular blocks - with courses'})];
        elseif C(j,2)==13
            tipology=[tipology; cell2table({'C1.1 Square cut stone wall - without courses'})];
        elseif C(j,2)==14
            tipology=[tipology; cell2table({'C1.2 Square cut stone wall - with courses'})];
        elseif C(j,2)==15
            tipology=[tipology; cell2table({'C2.0 Bricks'})];
        elseif C(j,2)==16
            tipology=[tipology; cell2table({'C2.1 Concrete blocks'})];
        end 
    end    
elseif agg==1
    % With aggregation
    for j=1:size(C,1)
        if C(j,2)==1
            tipology=[tipology; cell2table({'Irregular masonry - Without courses'})];
        elseif C(j,2)==2
            tipology=[tipology; cell2table({'Irregular masonry - With courses'})];
        elseif C(j,2)==3
            tipology=[tipology; cell2table({'Rough-hewn masonry - Without courses'})];
        elseif C(j,2)==4
            tipology=[tipology; cell2table({'Rough-hewn masonry - With courses'})];
        elseif C(j,2)==5
            tipology=[tipology; cell2table({'Square cut stone - Without courses'})];
        elseif C(j,2)==6
            tipology=[tipology; cell2table({'Square cut stone - With courses'})];
        elseif C(j,2)==7
            tipology=[tipology; cell2table({'Blocks'})];
        end 
    end
end
Output=[Output tipology];
Output.Properties.VariableNames={'No. of Builidngs','Main Resisting System'};

w=1;
% PARAMETER 1: Number of stories
if nnz(PAR(:)==1)>0
    par1={'1 Storey','2 Stories','3 Stories','4 Stories','5 Stories','6 Stories','7 Stories','8 Stories','9 Stories','10 Stories','11 Stories','12 Stories'};
    stories=[];
    for j=1:size(C,1)
        stories=[stories; cell2table(par1(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    stories.Properties.VariableNames={'Number of stories'};
    Output=[Output stories];
    w=w+2;
end

% PARAMETER 2: Construction period
if nnz(PAR(:)==2)>0
    if agg==0
        % Without Parameters Aggregation
        par2={'< 1860','1861-1919','1919-1945','1946-1961','1962-1971','1972-1975','1976-1981','1982-1986','1987-1991','1992-1996','1997-2001','2002-2008','2009-2011','> 2011'};
    elseif agg==1
        % With Parameters Aggregation
        par2={'< 1860','1861-1919','1919-1945','1945-1981','Post 1981'};
    end
    year=[];
    for j=1:size(C,1)
        year=[year; cell2table(par2(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    year.Properties.VariableNames={'Construction period'};
    Output=[Output year];
    w=w+2;
end

% PARAMETER 3: RING BEAMS AND TIE RODS
if nnz(PAR(:)==3)>0
    par3={'Present','Absent'};
    ringtie=[];
    for j=1:size(C,1)
        ringtie=[ringtie; cell2table(par3(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    ringtie.Properties.VariableNames={'Ring beams and Tie rods'};
    Output=[Output ringtie];
    w=w+2;
end

% PARAMETER 4: SLAB TYPE
if nnz(PAR(:)==4)>0
    if agg==0
        % Without Parameters Aggregation
        par4={'Timber beam slabs with terracotta tiles','Timber beam slabs with single wooden board','Steel beam slabs with brick arches', ...
            'Timber beam slabs with double wooden boards','SAP slabs','Steel beam slabs with flat brick tiles', ...
            'RC slabs','Slabs with inverted T beams prefabricated','Slabs with inverted T beams cast in situ'};
    elseif agg==1
        % With Parameters Aggregation
        par4={'Flexible','Semi-rigid','Rigid'};
    end
    slab=[];
    for j=1:size(C,1)
        slab=[slab; cell2table(par4(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    slab.Properties.VariableNames={'Slab type'};
    Output=[Output slab];
    w=w+2;
end

% PARAMETER 5: PRESENCE OF MIXED STRUCTURES
if nnz(PAR(:)==5)>0
    par5={'RC (or other framed structures) upwards extension on masonry','Masonry upwards extension on RC (or other framed structures)', ...
        'In-plan RC extension','Perimeter masonry and internal RC columns','Perimeter masonry an external columns','Confined masonry'};
    mixstruct=[];
    for j=1:size(C,1)
        mixstruct=[mixstruct; cell2table(par5(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    mixstruct.Properties.VariableNames={'Presence of mixed structure'};
    Output=[Output mixstruct];
    w=w+2;
end

% PARAMETER 6: REGULARITY IN PLAN
if nnz(PAR(:)==6)>0
    par6={'Regular','Moderately regular','Irregular'};
    regplan=[];
    for j=1:size(C,1)
        regplan=[regplan; cell2table(par6(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    regplan.Properties.VariableNames={'Regularity in plan'};
    Output=[Output regplan];
    w=w+2;
end

% PARAMETER 7: REGULARITY IN ELEVATION
if nnz(PAR(:)==7)>0
    par7={'Regular','Moderately regular','Irregular'};
    regelev=[];
    for j=1:size(C,1)
        regelev=[regelev; cell2table(par7(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    regelev.Properties.VariableNames={'Regularity in elevation'};
    Output=[Output regelev];
    w=w+2;    
end

% PARAMETER 8: ROOF MATERIAL
if nnz(PAR(:)==8)>0
    if agg==0
        % Without Parameters Aggregation 
        par8={'Wood','Steel','RC','Masonry'};        
    elseif agg==1
        % With Parameters Aggregation
        par8={'Light','Heavy'};
    end
    roof=[];
    for j=1:size(C,1)
        roof=[roof; cell2table(par8(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    roof.Properties.VariableNames={'Roof material'};
    Output=[Output roof];
    w=w+2;    
end

% PARAMETER 9: ROOF SPREAD
if nnz(PAR(:)==9)>0
    par9={'Yes','No'};
    spread=[];
    for j=1:size(C,1)
        spread=[spread; cell2table(par9(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    spread.Properties.VariableNames={'Roof spread'};
    Output=[Output spread];
    w=w+2;    
end

% PARAMETER 10: STRUCTURAL STRENGTHENING
if nnz(PAR(:)==10)>0
    par10={'Local strengthening','Seismic upgrade','Retrofit','None'};
    strength=[];
    for j=1:size(C,1)
        strength=[strength; cell2table(par10(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    strength.Properties.VariableNames={'Structural strengthening'};
    Output=[Output strength];
    w=w+2;
end

% PARAMETER 11: VAULT TYPE
if nnz(PAR(:)==11)>0
    par11={'Barrel vault','Barrel vault with bezels','Barrel vault with pavillon heads','Keel vault','Cloister vault','Groin vault','Bohemian vault','Fan vault'};
    vault=[];
    for j=1:size(C,1)
        vault=[vault; cell2table(par11(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    vault.Properties.VariableNames={'Vault type'};
    Output=[Output vault];
    w=w+2;    
end

% PARAMETER 12: AVERAGE STOREY HEIGHT
if nnz(PAR(:)==12)>0
    par12={'< 2.5 m','2.5-3.5 m','3.5-5 m','> 5 m'};
    avgheight=[];
    for j=1:size(C,1)
        avgheight=[avgheight; cell2table(par12(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    avgheight.Properties.VariableNames={'Average storey height'};
    Output=[Output avgheight];
    w=w+2;     
end

% PARAMETER 13: AVERAGE GROUND FLOOR HEIGHT
if nnz(PAR(:)==13)>0
    par13={'< 2.5 m','2.5-3.5 m','3.5-5 m','> 5 m'};
    avg1height=[];
    for j=1:size(C,1)
        avg1height=[avg1height; cell2table(par13(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    avg1height.Properties.VariableNames={'Average ground floor height'};
    Output=[Output avg1height];
    w=w+2; 
end

% PARAMETER 14: PERCENTAGE OF OPENINGS
if nnz(PAR(:)==14)>0
    par14={'< 10%','10-19%','20-29%','30-50%','>50%'};
    openings=[];
    for j=1:size(C,1)
        openings=[openings; cell2table(par14(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    openings.Properties.VariableNames={'Percentage of openings'};
    Output=[Output openings];
    w=w+2;     
end

% PARAMETER 15: PERCENTAGE OF OPENINGS AT THE GROUND FLOOR
if nnz(PAR(:)==15)>0
    par15={'< 10%','10-19%','20-29%','30-50%','>50%'};
    openings1=[];
    for j=1:size(C,1)
        openings1=[openings1; cell2table(par15(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    openings1.Properties.VariableNames={'Percentage of openings at the ground floor'};
    Output=[Output openings1];
    w=w+2; 
end

% PARAMETER 16: AVERAGE FLOOR AREA
if nnz(PAR(:)==16)>0
    par16={'50 m^2','70 m^2','100 m^2','130 m^2','170 m^2','230 m^2','300 m^2','400 m^2','500 m^2','650 m^2','900 m^2','1200 m^2','1600 m^2','2200 m^2','3000 m^2','> 3000 m^2'};
    area=[];
    for j=1:size(C,1)
        area=[area; cell2table(par16(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    area.Properties.VariableNames={'Average floor area'};
    Output=[Output area];
    w=w+2;     
end

% PARAMETER 17: AVERAGE WALL THICKNESS AT THE GROUND FLOOR
if nnz(PAR(:)==17)>0
    if agg==0
        par17={'< 40 cm','40-60 cm','> 60 cm'};
    elseif agg==1
        par17={'< 60 cm','> 60 cm'};
    end
    thick=[];
    for j=1:size(C,1)
        thick=[thick; cell2table(par17(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    thick.Properties.VariableNames={'Average wall thickness at the ground floor'};
    Output=[Output thick];
    w=w+2; 
end

% PARAMETER 18: AVERAGE DISTANCE BETWEEN WALLS
if nnz(PAR(:)==18)>0
    if agg==0
        par18={'< 4 m','4-5 m','5-6 m','> 6 m'};
    elseif agg==1
        par18={'< 5 m','> 5 m'};
    end
    distance=[];
    for j=1:size(C,1)
        distance=[distance; cell2table(par18(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end    
    distance.Properties.VariableNames={'Average distance between walls'};
    Output=[Output distance];
    w=w+2; 
end

% PERCENTAGES
perc=[table(num2str(C(:,end-1),'%.2f'),'VariableNames',{'Relevance (%)'}) table(num2str(C(:,end),'%.2f'),'VariableNames',{'Cumulative Relevance (%)'})];
Output=[Output perc];

ass=input('Perform sub-typologies assimilation? (Yes/No)\n','s');
if strcmpi(ass,'Yes')
    % Performing sub-typologies assimilation
    clearvars -except Error Output PAR PAR_selected C TOT
    openvar('Output')
    fprintf('\nProcess is paused. Please click "Continue" or type "dbcont" to resume.\n\n')
    keyboard
    
    % MATRIX C
    % Column 1: no. Buildings in Compartment
    % Column 2: ID Typology
    % Columns >2: Selected Parameters
    % Column end-1: Percentage on dataset
    % Column end: Cumulative percentage on dataset
    
    ril=input('Please provide number of chosen sub-typologies.\n');
    
    loc_tresh=input('Please provide local treshold for assimilation.\n');
    % How much the value for each parameter can vary from the one of the relevant sub-typology
    glob_tresh=input('Please provide global treshold for assimilation.\n');
    if glob_tresh<loc_tresh
        while glob_tresh<loc_tresh
            glob_tresh=input('Global treshold must be higher or equal than local treshold.\nPlease indicate global treshold for assimilation.\n');
        end
    end
    
    Graph=zeros(glob_tresh+1,2*ril);
    % Matrix to plot graphs of sub-typologies assimilation
    % Two columns per buildings selected as relevant
    % Odd columns, first row = Number of buildings for relevant sub-typology
    Graph(1,1:2:end)=C(1:ril,1)';
    % Even columns, second to last row = Number of buildings for assimilated sub-typology

    Assimilated=[];
    % Matrix of assimilated sub-typologies
    Discarded=[];
    % Matrix of NON-assimilated sub-typologies
    
    New_Perc=zeros(ril,2);
    % New relevance of sub-typology after assimilation    
    
    for j=1:ril
        % for each relevant sub-typology
        relevant=C(j,2);
        for w=1:length(PAR_selected)/2
            relevant=[relevant find(C(j,PAR_selected(2*w-1):PAR_selected(2*w)))];
        end
        
        for k=ril+1:size(C,1)
            % for each discarded sub-typology
            discarded=C(k,2);
            for w=1:length(PAR_selected)/2
                discarded=[discarded find(C(k,PAR_selected(2*w-1):PAR_selected(2*w)),1)];
            end
            
            diff=abs(relevant-discarded);
            % Row of differences
            % 1st column = differences in masonry type
            % 2nd to last column = differences in parameters 
            
            if diff(1)==0 && sum(diff(2:end))<=glob_tresh && max(diff(2:end))<=loc_tresh
                % Sub-typology can be assimilated
                Assimilated=[Assimilated;array2table(k) Output(k,:) array2table(j)];
                % Adding the number of relevant sub-typology it is assimilated to
                New_Perc(j,1)=New_Perc(j,1)+C(k,end-1);
                % Keeping track of sub-typology relevance after assimilation
                Graph(sum(diff(2:end))+1,2*j)=Graph(sum(diff(2:end))+1,2*j)+C(k,1);
                % Adding number of assimilated buildings for Graph
                C(k,:)=999;
                % Assimilated sub-typology is "removed" to prevent it from being assimilated again                
            end
        end
    end
    
    New_Perc(:,1)=New_Perc(:,1)+C(1:ril,end-1);
    % Add assimilated to original relevance
    New_Perc(:,2)=cumsum(New_Perc(:,1));
    Relevant=[Output(1:ril,:) array2table(New_Perc)];
    % Matrix of selected sub-typologies
    Relevant.Properties.VariableNames{'Relevance (%)'}='Original Relevance (%)';
    Relevant.Properties.VariableNames{'New_Perc1'}='Assimilated Relevance (%)';
    Relevant.Properties.VariableNames{'New_Perc2'}='Cumulative Assimilated Relevance (%)';
    
    % Matrix of assimilated sub-typologies
    Assimilated.Properties.VariableNames{'k'}='#';
    Assimilated.Properties.VariableNames{'j'}='Sub-typology Assimilated To';
    Assimilated=removevars(Assimilated,{'Cumulative Relevance (%)'});
    Assimilated=[Assimilated array2table(cumsum(str2double(table2cell(Assimilated(:,end-1)))))];
    Assimilated=movevars(Assimilated,'Var1','After','Relevance (%)');
    Assimilated.Properties.VariableNames{'Var1'}='Cumulative Relevance (%)';
    
    % Matrix of discarded sub-typologies
    disc=(ril+1:size(C,1))';
    [~,match]=ismember(table2array(Assimilated(:,1)),disc);
    disc(match,:)=[];
    Discarded=[array2table(disc) Output(disc,:)];
    Discarded.Properties.VariableNames{'disc'}='#';
    Discarded=removevars(Discarded,{'Cumulative Relevance (%)'});
    Discarded=[Discarded array2table(cumsum(str2double(table2cell(Discarded(:,end)))))];
    Discarded.Properties.VariableNames{'Var1'}='Cumulative Relevance (%)';    
    
    % Plotting Graphs
    n_row=2;
    % Number of rows in Figure
    n_col=3;
    % Number of columns in Figure
    w=1;
    k=1;
    j=1;
    while j<=ril 
        for k=1:n_row*n_col
            figure(w)
            if j>ril
                break
            end
            subplot(n_row,n_col,k)
            bar(Graph(:,2*j-1:2*j),'stacked'), grid on
            title(sprintf('Sub-typology %d',j))
            subtitle([sprintf('Frequency %.2f%%',Graph(1,2*j-1)/TOT*100),sprintf(' (Assimilated %.2f%%)',sum(sum(Graph(:,2*j-1:2*j)))/TOT*100)])
            yticks([0:sum(sum(Graph(:,2*j-1:2*j)))/10:sum(sum(Graph(:,2*j-1:2*j)))])
            ylim([0 sum(sum(Graph(:,2*j-1:2*j)))])
            yticklabels(num2cell(0:10:100))
            xticklabels(num2cell(0:1:glob_tresh))
            legend('Original','Assimilated','Location','northeast')
            xlabel('Difference')
            ylabel('Percentage')
            set(gca,'FontName','Arial')
            j=j+1;
        end
        w=w+1;
    end
    
    openvar('Discarded')
    openvar('Assimilated')
    openvar('Relevant')
    fprintf('\nAssimilation is complete. Please check Graphs.\n')
    clearvars -except Error Output Discarded Relevant Assimilated Graph      
    
elseif strcmpi(ass,'No')
    % Assimilation is NOT performed
    fprintf('\nPlease check Output and Errors.\n')
    clearvars -except Error Output    
end

fprintf('\nEND OF PROCESS.\n')