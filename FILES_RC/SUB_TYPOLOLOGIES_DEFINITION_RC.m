%% AUTOMATED SUBTYPOLOGIES DEFINITION BASED ON CARTIS DATA
% Alberto Basaglia, University "G. D'Annunzio" of Chieti-Pescara (2021)

% REINFORCED CONCRETE 

clear all
close all
clc

%% SUPPORTING FILES

Municipalities=readtable('Municipalities.txt','PreserveVariableNames',true);
% List of Italian Municipalities (as of JAN 2021)
Provinces=readtable('Provinces.txt','PreserveVariableNames',true);
% List of Italian Provinces (as of MARCH 2021)
% Province is used to define the seismic classification in the selected periods
% Periods: Pre-WWII, Post-WWII, DM'72, DDM'81-'85, DM'92-96, NTC2008-2018
% Possible entries: Gravity OR Seismic
% Distinction has to be done MANUALLY (supporting sofwtare, ECS-it)

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

% Conversion RC typologies in numbers (table to double)
tipology=CARTIS(:,8);
tipology(:,1)=regexprep(tipology{:,1},'','');
tipology=(strcmp(tipology{:, 1},'A - Telai Tamponati')*1+strcmp(tipology{:,1},'B - Telai Con Travi')*2 ...
       +strcmp(tipology{:,1},'C - Telai Travi Spessore')*3+strcmp(tipology{:,1},'D - Telai Travi Su Perimetro')*4 ...
       +strcmp(tipology{:,1},'E - Telai Travi E NucleiCA')*5+strcmp(tipology{:,1},'F - Setti')*6 ...
       +strcmp(tipology{:,1},'G - Telai E Setti')*7);
% NOTE: typo in Typology E (from pgAdmin database)
A=[A tipology];

disp('Parameters for disaggregation (pgAdmin MARCH 2021).')
disp(' PAR 1: Number of stories')
disp(' PAR 2: Construction period')
disp(' PAR 3: Seismic joints')
disp(' PAR 4: Unidirectional frames')
disp(' PAR 5: Ground floor infill walls')
disp(' PAR 6: Soft storey upper floors')
disp(' PAR 7: Presence of SAP slabs (or similar)')
disp(' PAR 8: Regularity in plan')
disp(' PAR 9: Regularity in elevation')
disp(' PAR 10: Roof material')
disp(' PAR 11: Roof spread')
disp(' PAR 12: Structural strengthening')

PAR=input('\nInsert number of selected parameters in [], space delimited:\n');
clc

% Columns     1     2     3     4     5     6     7     8     9     10    11   12
PAR_columns=[9 20 22 35 37 38 40 41 43 45 47 48 50 51 53 55 57 59 61 64 66 67 69 72];
PAR_selected=[];

for j=1:12
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
% Fix Parameter 3: Seismic joints
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
% Fix Parameter 4: Unidirectional Frames
if nnz(PAR(:)==4)>0
    pos=find(PAR(:)==4);
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
% Fix Parameter 7: Presence of SAP slabs (or similar)
if nnz(PAR(:)==7)>0
    pos=find(PAR(:)==7);
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
% Fix Parameter 12: Structural strengthening
if nnz(PAR(:)==12)>0
    pos=find(PAR(:)==12);
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
            % CASE 3: Data is present and consistent: proceed with disaggregation
            B_dis=zeros(prod(z),size(B_temp,2));
            % Temporary matrix with disaggregated compartment
            B_dis(:,2)=(B_temp(k,3));
            % ID Province is copied (MOVED TO 2nd COLUMN)
            B_dis(:,3)=(B_temp(k,4));
            % ID Typology is copied
            znot=find(B_temp(k,:));
            % NOTE: 4th Column of B_dis (and B) will be used to define if design
            % followed Gravitational or Seismic loads IF Parameter 2 has
            % been selected for disaggregation
            
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
            comb=allcomb(z1,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12);
            comb(:,~any(comb,1))=[];
            % Columns with 0 are removed
        
            for l=1:prod(z)
                B_dis(l,comb(l,:))=1;
                % A "1" is put in designated columns (see comb) to mark the value
                % of selected parameters
                B_dis(l,1)=ceil(B_temp(k,1)*prod(B_temp(k,comb(l,:))./B_temp(k,1)));
                % No. of Buildings in disaggregated Compartment
                % (equally distributed according to distribution among selected parameter)
            end
        
            B=[B; B_dis];
            % Output    
            
        end
       
    end
    % Cycle repeats for next Municipality
end

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
    fprintf('The disaggregation is complete. %d buildings were considered (%.2f%% of the total).\n\n',sum(B(:,1)),min(100,(sum(B(:,1))/sum(table2array(CARTIS(:,3))))*100));
end

% MATRIX B
% Column 1: no. Buildings in Compartment
% Column 2: ID Province
% Column 3: ID Typology
% Column 4: Gravity/Seismic loads (IF Parameter 2 has been selected)
% Columns >4: Selected Parameters

% ASSOCIATING PROVINCE AND CONSTRUCTION PERIOD TO DESIGN LOAD
% B, 4th Column: 1 = Gravity loads
% B, 4th Column: 2 = Seismic loads
if nnz(PAR==2)>0
    pos=PAR_selected(2*find(PAR==2)-1);
    % Find position of first column of Parameter 2 in B
    for j=1:size(B,1)               
        % 6 Intervals: Pre-1945
        %              Post-1945
        %              DM'72
        %              DM'80 and '85
        %              DM'82 and '96
        %              Post-NTC2008
        period=find([nnz(B(j,pos:pos+2)) nnz(B(j,pos+3:pos+4)) nnz(B(j,pos+5:pos+6)) nnz(B(j,pos+7:pos+8)) nnz(B(j,pos+9:pos+11)) nnz(B(j,pos+12:pos+13))]);
        % Finding in which interval the Buildings in the j-th row were made
        design=strcmp(Provinces{B(j,2),3:8},'Seismic');
        design=design(period);
        if design==1
            % Seismic Loads
            B(j,4)=2;
        else
            % Gravity Loads
            B(j,4)=1;
        end
    end
end

%% PARAMETERS' AGGREGATION (OPTIONAL)

agg=input('Perform parameters aggregation? (Yes/No)\n','s');
if strcmpi(agg,'Yes')
    agg=1;
elseif strcmpi(agg,'No')
    agg=0;
end

% Performing aggregation
if agg==1        
    % Parameter 1: Number of stories
    if nnz(PAR==1)>0
        pos=PAR_selected(2*find(PAR==1)-1);
        % Find position of first column of Parameter 1 in B
        B=[B(:,1:pos+7) sum(B(:,pos+8:pos+11),2) B(:,pos+12:end)];
        % Buildings with 9 or more stories -> Tall Buildings
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-3;
        % Updating of PAR_selected given the condensed number of columns
    end
    
    % Parameter 2: Construction period
    if nnz(PAR==2)>0
        pos=PAR_selected(2*find(PAR==2)-1);
        % Find position of first column of Parameter 2 in B
        B=[B(:,1:pos-1) sum(B(:,pos:pos+2),2) sum(B(:,pos+3:pos+4),2) sum(B(:,pos+5:pos+6),2) sum(B(:,pos+7:pos+8),2) sum(B(:,pos+9:pos+11),2) sum(B(:,pos+12:pos+13),2) B(:,pos+14:end)];
        % Construction periods -> 6 intervals defined in Province
        %                         (evolution of building codes) 
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-8;
        % Updating of PAR_selected given the condensed number of columns
    end
    
    % Parameter 5: Ground floor infill walls
    if nnz(PAR==5)>0
        pos=PAR_selected(2*find(PAR==5)-1);
        % Find position of first column of Parameter 5 in B
        B=[B(:,1:pos) sum(B(:,pos+1:pos+2),2) B(:,pos+3:end)];
        % Infills with irregular disposition or absent -> Infills absent
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-1;
        % Updating of PAR_selected given the condensed number of columns                
    end    
   
    % Parameter 8: Regularity in plan
    if nnz(PAR==8)>0
        pos=PAR_selected(2*find(PAR==8)-1);
        % Find position of first column of Parameter 8 in B
        B=[B(:,1:pos-1) sum(B(:,pos:pos+1),2) B(:,pos+2:end)];
        % Buildings regular or moderately regular in plan -> Regular in plan
        PAR_selected(PAR_selected>pos)=PAR_selected(PAR_selected>pos)-1;
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
    % ID Province and Typology are omitted (as the former is included in the
    % design type, 4th column, and in the cicle, i.e. j)
    while isempty(C1)==0
        C1diff=C1-[0 C1(1,2:end)];
        % To look for other rows with buildings having the same
        % characteristics of the first 
        % NOTE: first column is not considered, i.e. no. of buildings
        d1=find(all(C1diff(:,2:end)==0,2));
        % Rows of buildings with the same characteristics
        C=[C; [sum(C1diff(d1,1)) j C1(1,2:end)]];
        % ID Province is dropped
        C1(d1,:)=[];
        % Removing rows summed
    end
end
PAR_selected=PAR_selected-1;
% As ID Province has been dropped

% Sorting and adding percentages
C=[C (C(:,1)./sum(C(:,1)))*100];
C=sortrows(C,size(C,2),'descend');
C=[C cumsum(C(:,end))];

TOT=sum(C(:,1));

% MATRIX C
% Column 1: no. Buildings in Compartment
% Column 2: ID Typology
% Column 3: Gravity/Seismic loads (or 0 if Parameter 2 has NOT been selected)
% Columns >3: Selected Parameters
% Column end-1: Percentage on dataset
% Column end: Cumulative percentage on dataset

%% OUTPUT

Output=[];

% NO. OF BUILDINGS AND TYPOLOGY
Output=[Output array2table(C(:,1))];
Output.Properties.VariableNames={'No. of Builidngs'};
tipology=[];
for j=1:size(C,1)
    if C(j,2)==1
        tipology=[tipology; cell2table({'A. Frames with sturdy infill walls'})];
    elseif C(j,2)==2
        tipology=[tipology; cell2table({'B. Frames with exposed beams and non-sturdy infill walls'})];
    elseif C(j,2)==3
        tipology=[tipology; cell2table({'C. Frames with flat beams and non-sturdy or absent infill walls'})];
    elseif C(j,2)==4
        tipology=[tipology; cell2table({'D. Frames with exposed beam on the perimeter, non-sturdy or absent infill walls and inner flat beams'})];
    elseif C(j,2)==5
        tipology=[tipology; cell2table({'E. Dual system, frames with exposed beams and inner RC shear walls '})];
    elseif C(j,2)==6
        tipology=[tipology; cell2table({'F. Shear walls'})];
    elseif C(j,2)==7
        tipology=[tipology; cell2table({'G. Dual system, frames with flat beams and inner shear walls'})];
    end 
end
Output=[Output tipology];
Output.Properties.VariableNames={'No. of Builidngs','Main Resisting System'};

w=1;
% PARAMETER 1: Number of stories
if nnz(PAR(:)==1)>0
    if agg==0
        % Without Parameters Aggregation
        par1={'1 Storey','2 Stories','3 Stories','4 Stories','5 Stories','6 Stories','7 Stories','8 Stories','9 Stories','10 Stories','11 Stories','12 Stories'};
    elseif agg==1
        % With Parameters Aggregation
        par1={'1 Storey','2 Stories','3 Stories','4 Stories','5 Stories','6 Stories','7 Stories','8 Stories','>8 Stories'};
    end
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
       par2={'Before World War II','After World War II','DM 72','DM 81 or 85','DM 92 or 96','Post NTC 2008'};
   end
    year=[];
    design=[];
    for j=1:size(C,1)
        year=[year; cell2table(par2(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
        if C(j,3)==1
            design=[design; cell2table({'Gravity'})];
        else
            design=[design; cell2table({'Seismic'})];
        end
    end
    year.Properties.VariableNames={'Construction period'};
    design.Properties.VariableNames={'Design load'};
    Output=[Output year design];
    w=w+2;
end

% PARAMETER 3: Seismic joints
if nnz(PAR(:)==3)>0
    par3={'Code compliant','Non-compliant'};
    joints=[];
    for j=1:size(C,1)
        joints=[joints; cell2table(par3(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    joints.Properties.VariableNames={'Seismic joints'};
    Output=[Output joints];
    w=w+2;
end

% PARAMETER 4: Unidirectional frames
if nnz(PAR(:)==4)>0
    par4={'Yes','No'};
    frames=[];
    for j=1:size(C,1)
        frames=[frames; cell2table(par4(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    frames.Properties.VariableNames={'Unidirectional frames'};
    Output=[Output frames];
    w=w+2;
end

% PARAMETER 5: Ground floor infill walls
if nnz(PAR(:)==5)>0
    if agg==0
        % Without Parameters Aggregation
        par5={'Regular configuration','Irregular configuration','Absent'};
    elseif agg==1
        % With Parameters Aggregation
        par5={'Regular configuration','Absent'};
    end
    infills=[];
    for j=1:size(C,1)
        infills=[infills; cell2table(par5(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    infills.Properties.VariableNames={'Ground floor infill walls'};
    Output=[Output infills];
    w=w+2;
end

% PARAMETER 6: Soft storey upper floors
if nnz(PAR(:)==6)>0
    par6={'Yes','No'};
    softstorey=[];
    for j=1:size(C,1)
        softstorey=[softstorey; cell2table(par6(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    softstorey.Properties.VariableNames={'Soft storey upper floors'};
    Output=[Output softstorey];
    w=w+2;
end

% PARAMETER 7: Presence of SAP slabs (or similar)
if nnz(PAR(:)==7)>0
    par7={'Yes','No'};
    sapslabs=[];
    for j=1:size(C,1)
        sapslabs=[sapslabs; cell2table(par7(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    sapslabs.Properties.VariableNames={'Presence of SAP slabs (or similar)'};
    Output=[Output sapslabs];
    w=w+2;
end

% PARAMETER 8: Regularity in plan
if nnz(PAR(:)==8)>0
   if agg==0
        % Without Parameters Aggregation
        par8={'Regular','Moderately regular','Irregular'};
   elseif agg==1
        % With Parameters Aggregation
       par8={'Regular','Irregular'};
   end
    regplan=[];
    for j=1:size(C,1)
        regplan=[regplan; cell2table(par8(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    regplan.Properties.VariableNames={'Regularity in plan'};
    Output=[Output regplan];
    w=w+2;
end    

% PARAMETER: Regularity in elevation
if nnz(PAR(:)==9)>0
    par9={'Regular','Moderately regular','Irregular'};
    regelev=[];
    for j=1:size(C,1)
        regelev=[regelev; cell2table(par9(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    regelev.Properties.VariableNames={'Regularity in elevation'};
    Output=[Output regelev];
    w=w+2;
end

% PARAMETER 10: Roof material
if nnz(PAR(:)==10)>0
    par10={'Wood','Steel','RC','Masonry'};
    roof=[];
    for j=1:size(C,1)
        roof=[roof; cell2table(par10(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    roof.Properties.VariableNames={'Roof material'};
    Output=[Output roof];
    w=w+2;
end
    
% PARAMETER 11: Roof spread
if nnz(PAR(:)==11)>0
    par11={'Yes','No'};
    spread=[];
    for j=1:size(C,1)
        spread=[spread; cell2table(par11(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    spread.Properties.VariableNames={'Roof spread'};
    Output=[Output spread];
    w=w+2;
end

% PARAMETER 12: Structural strengthening
if nnz(PAR(:)==12)>0
    par12={'Local strengthening','Seismic upgrade','Retrofit','None'};
    strength=[];
    for j=1:size(C,1)
        strength=[strength; cell2table(par12(find(C(j,PAR_selected(w):PAR_selected(w+1)))))];
    end
    strength.Properties.VariableNames={'Structural strengthening'};
    Output=[Output strength];
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
    % Column 3: Gravity/Seismic loads (or 0 if Parameter 2 has NOT been selected)
    % Columns >3: Selected Parameters
    % Column end-1: Percentage on dataset
    % Column end: Cumulative percentage on dataset
    
    ril=input('Please select number of chosen sub-typologies.\n');
    
    loc_tresh=input('Please indicate local treshold for assimilation.\n');
    % How much the value for each parameter can vary from the one of the relevant sub-typology
    glob_tresh=input('Please indicate global treshold for assimilation.\n');
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
        if nnz(PAR==2)>0
            relevant=[C(j,2) C(j,3)];
        else
            relevant=C(j,2);
        end
        
        for w=1:length(PAR_selected)/2
            relevant=[relevant find(C(j,PAR_selected(2*w-1):PAR_selected(2*w)))];
        end
        
        for k=ril+1:size(C,1)
            % for each discarded sub-typology
            if nnz(PAR==2)>0
                discarded=[C(k,2) C(k,3)];
            else
                discarded=C(k,2);
            end
            
            for w=1:length(PAR_selected)/2
                discarded=[discarded find(C(k,PAR_selected(2*w-1):PAR_selected(2*w)),1)];
            end
            
            diff=abs(relevant-discarded);
            % Row of differences
            % 1st column = differences in main resisting system
            % 2nd column = differences in design type (ONLY IF CONSTRUCTION PERIOD IS SELECTED)
            % 3rd to last column = differences in parameters
            
            if nnz(PAR==2)>0 && diff(1)==0 && diff(2)==0 && sum(diff(3:end))<=glob_tresh && max(diff(3:end))<=loc_tresh
                % CASE 1: Parameter "Costruction Period" has been selected
                % Sub-typology can be assimilated
                Assimilated=[Assimilated;array2table(k) Output(k,:) array2table(j)];
                % Adding the number of relevant sub-typology it is assimilated to
                New_Perc(j,1)=New_Perc(j,1)+C(k,end-1);
                % Keeping track of sub-typology relevance after assimilation
                Graph(sum(diff(3:end))+1,2*j)=Graph(sum(diff(3:end))+1,2*j)+C(k,1);
                % Adding number of assimilated buildings for Graph
                C(k,:)=1;
                % Assimilated sub-typology is "removed" to prevent it from being assimilated again
            elseif nnz(PAR==2)==0 && diff(1)==0 && sum(diff(2:end))<=glob_tresh && max(diff(2:end))<=loc_tresh
                % CASE 2: Parameter "Costruction Period" has NOT been selected
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
            yticks([0:sum(sum(Graph(:,2*j-1:2*j)))/5:sum(sum(Graph(:,2*j-1:2*j)))])
            ylim([0 sum(sum(Graph(:,2*j-1:2*j)))])
            yticklabels(num2cell(0:20:100))
            xticklabels(num2cell(0:1:glob_tresh))
            legend('Original','Assimilated','Location','northeast')
            xlabel('Difference')
            ylabel('Percentage')
            set(gca,'FontName','Arial','FontSize',12)
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