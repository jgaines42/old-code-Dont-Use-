%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [are_hitting] = move_check_2(allDipeptide1,allDipeptide2,check_all, rest_backbone)
%
% Determines if two residues can hit each other
%
% Input: 
% allDipeptide1
% allDipeptide2
% check_all
% rest_backbone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [are_hitting] = move_check_2(allDipeptide1,allDipeptide2,check_all, rest_backbone)

overlap_amount = 0;
Back_Pos = cell2mat(rest_backbone(:,8:10));
%Extract Residue names
resiName1 = allDipeptide1{4,4};
resiName1(2:3) = lower(resiName1(2:3));
resiID1 = cell2mat(allDipeptide1(4,6));
resiName2 = allDipeptide2{4,4};
resiName2(2:3) = lower(resiName2(2:3));

%Val is always first and Met is always last
%{
if resiName2 == 'Val'
    temp = allDipeptide2;
    allDipeptide2 = allDipeptide1;
    allDipeptide1 = temp;
    resiName1 = allDipeptide1{4,4};
    resiName1(2:3) = lower(resiName1(2:3));
    resiID1 = cell2mat(allDipeptide1(4,6));
    resiName2 = allDipeptide2{4,4};
    resiName2(2:3) = lower(resiName2(2:3));
elseif resiName1 == 'Met'
    temp = allDipeptide2;
    allDipeptide2 = allDipeptide1;
    allDipeptide1 = temp;
    resiName1 = allDipeptide1{4,4};
    resiName1(2:3) = lower(resiName1(2:3));
    resiID1 = cell2mat(allDipeptide1(4,6));
    resiName2 = allDipeptide2{4,4};
    resiName2(2:3) = lower(resiName2(2:3));
end
%}
%Step size will be 10 degrees
max_Chi2_1 = 36;
max_Chi2_2 = 36;
if resiName1 == 'Phe'
    max_Chi2_1 = 36/2;
end
if resiName2 == 'Phe'
    max_Chi2_2 = 36/2;
end

resiID2 = cell2mat(allDipeptide2(4,6));
Position_1 = cell2mat(allDipeptide1(:,8:10));
Position_2 = cell2mat(allDipeptide2(:,8:10));
Xtal1 = Position_1;
Xtal2 = Position_2;

%Set up Residue 1
res_name = resiName1;
switch_residue_setup;
iChi1Array1 = iChi1Array;
moveAtomID2_1=moveAtomID2; %%%%%%%% rotate chi2
InitChi1_1=calcDA2(iChi1Array,Position_1); %%%%%%%%
InitChi1_1=mod(real(InitChi1_1),360);
if DOF > 1
    iChi2Array1 = iChi2Array;
    moveAtomID_1=moveAtomID; %%%%%%%% rotate chi1
    InitChi2_1=calcDA2(iChi2Array,Position_1); %%%%%%%%
    InitChi2_1=mod(real(InitChi2_1),360);
end
if DOF > 2
    iChi3Array1 = iChi3Array;
    moveAtomID3_1 = moveAtomID3;
    InitChi3_1=calcDA2(iChi3Array,Position_1); %%%%%%%%
    InitChi3_1=mod(real(InitChi3_1),360);
    
end
numAtom1 =numAtom;
DOF1 =DOF;
all_move_1 = moveAtomID2; %Which atoms wil lmove

%% Set up Residue 2
res_name = resiName2;
switch_residue_setup;
iChi1Array2 = iChi1Array;
moveAtomID2_2=moveAtomID2; %%%%%%%% rotate chi1
InitChi1_2=calcDA2(iChi1Array2,Position_2); %%%%%%%%
InitChi1_2=mod(real(InitChi1_2),360);
if DOF > 1
    iChi2Array2 = iChi2Array;
    moveAtomID_2=moveAtomID; %%%%%%%% rotate chi2
    InitChi2_2=calcDA2(iChi2Array2,Position_2); %%%%%%%%
    InitChi2_2=mod(real(InitChi2_2),360);
end
if DOF > 2
    iChi3Array2 = iChi3Array;
    moveAtomID3_2 = moveAtomID3;
    InitChi3_2=calcDA2(iChi3Array2,Position_2); %%%%%%%%
    InitChi3_2=mod(real(InitChi3_2),360);
end
numAtom2 =numAtom;
DOF2 =DOF;
all_move_2 = moveAtomID2;

%Check that you have the right number of atoms
if size(Position_1,1) ~= numAtom1 || size(Position_2,1) ~= numAtom2
    are_hitting = -1;
else
    %
    %Make clash list
    res_clash= [];
    for clash_res = 1:size(moveAtomID2_2,2)
        res_clash = [res_clash;moveAtomID2_1', repmat(moveAtomID2_2(clash_res),size(moveAtomID2_1,2),1)];
    end
    
    res_clash(:,3) = cell2mat(allDipeptide1(res_clash(:,1),12))  + cell2mat(allDipeptide2(res_clash(:,2),12));
    res_clash(:,4) = res_clash(:,3).^2;
    

    
    %% If both are Val
    subtract_array_1_1 = repmat(Position_1(iChi1Array1(2),:),numAtom1,1);
    delta_term_1_1 =  pi*sign(InitChi1_1)*InitChi1_1/180;
    subtract_array_1_2 = repmat(Position_2(iChi1Array2(2),:),numAtom2,1);
    delta_term_1_2 =  pi*sign(InitChi1_2)*InitChi1_2/180;
    
    if check_all == 1
        %Make clash list for backbone
        back_clash1 = zeros(size(moveAtomID2_1,2)*size(rest_backbone,1),4);
        for clash_res = 1:size(moveAtomID2_1,2)
            back_clash1((clash_res-1)*size(rest_backbone,1)+1 : clash_res*size(rest_backbone,1),1:2) = [repmat(moveAtomID2_1(clash_res),size(rest_backbone,1),1),[1:size(rest_backbone,1)]'];
        end
        back_clash1(:,3) = cell2mat(allDipeptide1(back_clash1(:,1),12))+cell2mat(rest_backbone(back_clash1(:,2),12));
        back_clash1(:,4) = back_clash1(:,3).^2;
        
        dist = Position_1(back_clash1(:,1),1:3)-Back_Pos(back_clash1(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - back_clash1(:,4);
        ind0 = clash < 100;
        back_clash1 = back_clash1(ind0,:);
        
        back_clash2 = zeros(size(moveAtomID2_2,2)*size(rest_backbone,1),4);
        for clash_res = 1:size(moveAtomID2_2,2)
            back_clash2((clash_res-1)*size(rest_backbone,1)+1 : clash_res*size(rest_backbone,1),1:2) = [repmat(moveAtomID2_2(clash_res),size(rest_backbone,1),1),[1:size(rest_backbone,1)]'];
        end
        back_clash2(:,3) = cell2mat(allDipeptide2(back_clash2(:,1),12))+cell2mat(rest_backbone(back_clash2(:,2),12));
        back_clash2(:,4) = back_clash2(:,3).^2;
        dist = Position_2(back_clash2(:,1),1:3)-Back_Pos(back_clash2(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - back_clash2(:,4);
        ind0 = clash < 100;
        back_clash2 = back_clash2(ind0,:);
        
        
        %% Do single residue of each, get sublist of things to check
        D1 = [];
         for chi1_1 = 1:36 % Rotate Chi1 of Residue 1
            Position_1 = Xtal1;
            setChi1 = chi1_1*10;
            Position_1 = Rotate_DA(Position_1, setChi1, subtract_array_1_1, delta_term_1_1, iChi1Array1, moveAtomID2_1);
            if DOF1 == 1
                dist = Position_1(back_clash1(:,1),1:3)-Back_Pos(back_clash1(:,2),1:3);
                distemp = sqrt(sum(dist.^2,2));
                clash = distemp - back_clash1(:,3);
                if min(clash)>=overlap_amount
                    D1 = [D1;setChi1];
                end
            elseif DOF1 >= 2
                Pos_b4_Chi2_1 = Position_1;
                subtract_array_2_1 = repmat(Position_1(iChi2Array1(2),:),numAtom1,1);
                delta_term_2_1 =  pi*sign(InitChi2_1)*InitChi2_1/180;
                
                for chi2_1 = 1:max_Chi2_1
                    Position_1 = Pos_b4_Chi2_1;
                    setChi2 = chi2_1*10;
                    Position_1 = Rotate_DA(Position_1, setChi2, subtract_array_2_1, delta_term_2_1, iChi2Array1, moveAtomID_1);
                    if DOF1 == 2
                        %Check for backbone clash of residue 1
                        dist = Position_1(back_clash1(:,1),1:3)-Back_Pos(back_clash1(:,2),1:3);
                        distemp = sqrt(sum(dist.^2,2));
                        clash = distemp - back_clash1(:,3);
                        if min(clash)>=overlap_amount
                            D1 = [D1;setChi1, setChi2];
                        end
                    else
                        Pos_b4_Chi3 = Position_1;
                         subtract_array_3_1 = repmat(Position_1(iChi3Array1(2),:),numAtom1,1);
                         delta_term_3_1 =  pi*sign(InitChi3_1)*InitChi3_1/180;
                         
                         for chi3_1 = 1:36 %Rotate Chi3 of Residue 2
                             Position_1 = Pos_b4_Chi3;
                             setChi3 = chi3_1*10;
                             Position_1 = Rotate_DA(Position_1, setChi3, subtract_array_3_1, delta_term_3_1, iChi3Array1, moveAtomID3_1);
                             
                             %Check for backbone clash of Residue 2
                             dist = Position_1(back_clash1(:,1),1:3)-Back_Pos(back_clash1(:,2),1:3);
                             distemp = sqrt(sum(dist.^2,2));
                             clash = distemp - back_clash1(:,3);
                             if min(clash)>=overlap_amount
                                 D1 = [D1;setChi1, setChi2, setChi3];
                             end
                         end
                    end
                        
                end
            end
         end
         
         D2 = [];
         for chi1_2 = 1:36 %Rotate Chi1 of Residue 2
             Position_2 = Xtal2;
             setChi1 = chi1_2*10;
             Position_2 = Rotate_DA(Position_2, setChi1, subtract_array_1_2, delta_term_1_2, iChi1Array2, moveAtomID2_2);
             
             if DOF2 == 1 % Residue 2 only has 1 DOF
                 % Check for backbone clash of residue 2
                 dist = Position_2(back_clash2(:,1),1:3)-Back_Pos(back_clash2(:,2),1:3);
                 distemp = sqrt(sum(dist.^2,2));
                 clash = distemp - back_clash2(:,3);
                 if min(clash)>=overlap_amount
                    D2 = [D2;setChi1];
                 end
             else
                 Pos_b4_Chi2 = Position_2;
                 subtract_array_2_2 = repmat(Position_2(iChi2Array2(2),:),numAtom2,1);
                 delta_term_2_2 =  pi*sign(InitChi2_2)*InitChi2_2/180;
                 
                 for chi2_2 = 1:max_Chi2_2 %Rotate Residue 2 Chi2
                     Position_2 = Pos_b4_Chi2;
                     setChi2 = chi2_2*10;
                     Position_2 = Rotate_DA(Position_2, setChi2, subtract_array_2_2, delta_term_2_2, iChi2Array2, moveAtomID_2);
                     if DOF2 == 2 %If Residue 2 only has 2 DOF
                         %Check for backbone clash of Residue 2
                         dist = Position_2(back_clash2(:,1),1:3)-Back_Pos(back_clash2(:,2),1:3);
                        distemp = sqrt(sum(dist.^2,2));
                        clash = distemp - back_clash2(:,3);
                         if min(clash)>=overlap_amount
                             D2 = [D2;setChi1, setChi2];
                         end
                     else %DOF == 3
                         Pos_b4_Chi3 = Position_2;
                         subtract_array_3_2 = repmat(Position_2(iChi3Array2(2),:),numAtom2,1);
                         delta_term_3_2 =  pi*sign(InitChi3_2)*InitChi3_2/180;
                         
                         for chi3_2 = 1:36 %Rotate Chi3 of Residue 2
                             Position_2 = Pos_b4_Chi3;
                             setChi3 = chi3_2*10;
                             Position_2 = Rotate_DA(Position_2, setChi3, subtract_array_3_2, delta_term_3_2, iChi3Array2, moveAtomID3_2);
                             
                             %Check for backbone clash of Residue 2
                             dist = Position_2(back_clash2(:,1),1:3)-Back_Pos(back_clash2(:,2),1:3);
                             distemp = sqrt(sum(dist.^2,2));
                             clash = distemp - back_clash2(:,3);
                             if min(clash)>=overlap_amount
                                 D2 = [D2;setChi1,setChi2,setChi3];
                             end
                         end
                     end
                 end
             end
         end

        %%
        are_hitting = 0;
        for loop_1 = 1:size(D1,1)
            Position_1 = Xtal1;
            setChi1 = D1(loop_1,1);
            Position_1 = Rotate_DA(Position_1, setChi1, subtract_array_1_1, delta_term_1_1, iChi1Array1, moveAtomID2_1);
            if DOF1 >1
                setChi2 = D1(loop_1,2);
                subtract_array_2_1 = repmat(Position_1(iChi2Array1(2),:),numAtom1,1);
                Position_1 = Rotate_DA(Position_1, setChi2, subtract_array_2_1, delta_term_2_1, iChi2Array1, moveAtomID_1);
                if DOF1 > 2
                    setChi3 = D1(loop_1,3);
                    subtract_array_3_1 = repmat(Position_1(iChi3Array1(2),:),numAtom1,1);
                    Position_1 = Rotate_DA(Position_1, setChi3, subtract_array_3_1, delta_term_3_1, iChi3Array1, moveAtomID3_1);
                end
            end
            
            for loop_2 = 1:size(D2,1)
                Position_2 = Xtal2;
                setChi1 = D2(loop_2,1);
                Position_2 = Rotate_DA(Position_2, setChi1, subtract_array_1_2, delta_term_1_2, iChi1Array2, moveAtomID2_2);
                if DOF2 >1
                    setChi2 = D2(loop_2,2);
                    
                    subtract_array_2_2 = repmat(Position_2(iChi2Array2(2),:),numAtom2,1);
                    delta_term_2_2 =  pi*sign(InitChi2_2)*InitChi2_2/180;
                    Position_2 = Rotate_DA(Position_2, setChi2, subtract_array_2_2, delta_term_2_2, iChi2Array2, moveAtomID_2);
                    if DOF2 > 2
                        setChi3 = D2(loop_2,3);
                        subtract_array_3_2 = repmat(Position_2(iChi3Array2(2),:),numAtom2,1);
                        delta_term_3_2 =  pi*sign(InitChi3_2)*InitChi3_2/180;
                        Position_2 = Rotate_DA(Position_2, setChi3, subtract_array_3_2, delta_term_3_2, iChi3Array2, moveAtomID3_2);
                    end
                end
                dist = Position_1(res_clash(:,1),1:3)-Position_2(res_clash(:,2),1:3);
                distemp = sum(dist.^2,2);
                clash = distemp - res_clash(:,4);
                ind_clash = find(clash < 0);
                if size(ind_clash,1) > 0
                    are_hitting = 1;
                    break;
                end
            end
            if are_hitting == 1;
                break
            end
        end
    else %just check initial
        dist = Position_1(res_clash(:,1),1:3)-Position_2(res_clash(:,2),1:3);
        distemp = sum(dist.^2,2);
        clash = distemp - res_clash(:,4);
        ind_clash = find(clash < 0);
        if size(ind_clash,1) > 0
            are_hitting = 1;
        else
            are_hitting = 0;
            
        end
        
        
    end
end
end
          
           