% ����stk����Starlink����
% ������ 2019.06.23
% function StarlinkAnalyze_COM
% ��STK�����½�����
% Requires STK, STK Coverage and STK Integration licenses.
clc
clear all
% % 
% % % Initialize 
disp('��STK���򣬴���һ������');
try
    % ��ȡ�����е�STKʵ�����
    uiapp = actxGetRunningServer('STK11.application');
    root = uiapp.Personality2;
    checkempty = root.Children.Count;
    if checkempty == 0
        % ���δ���ֳ������½�һ��
        uiapp.visible = 1;
        root.NewScenario('StarLink_COM');
        scenario = root.CurrentScenario;
    else
        % ������ִ��ŵĳ�����ѯ���Ƿ�ر��ؽ�
        rtn = questdlg({'Close the current scenario?',' ','(WARNING: If you have not saved your progress will be lost)'});
        if ~strcmp(rtn,'Yes')
            return
        else
            root.CurrentScenario.Unload
            uiapp.visible = 1;
            root.NewScenario('StarLink_COM');
            scenario = root.CurrentScenario;
        end
    end

catch
    % STKû���������½�ʵ����ȡ���
    uiapp = actxserver('STK11.application');
    root = uiapp.Personality2;
    uiapp.visible = 1;
    root.NewScenario('StarLink_COM');
    scenario = root.CurrentScenario;
end

% �趨����ʱ��
disp('�趨����ʱ��');
strBegTime = '27 May 2019 06:14:00.000';
strEndTime = '28 May 2019 06:14:00.000';
scenario.SetTimePeriod(strBegTime,strEndTime);
% �趨����������ʼʱ��
scenario.Epoch = strBegTime;
% �趨����ʱ��ص���ʼ�㣬�����õ���Connect����
root.ExecuteCommand('Animate * Reset');

% �½���������
disp('������������');
seedSat = scenario.Children.New('eSatellite', 'Sat');
seedSat.SetPropagatorType('ePropagatorJ4Perturbation');

% �������Ǹ���
% ע�����ﵥλ��mexConnect�в�ͬ���볤���ǹ���,�Ƕ��Ƕ�
a=6928.137;  e=0.0; i=53.0;w=0; Raan=160; M=0;
propagator = seedSat.Propagator; 
propagator.StartTime = scenario.StartTime;  % Set to scenario start time
propagator.StopTime = scenario.StopTime;    % Set to scenario stop time
propagator.Step = 60.0;
propagator.InitialState.Representation.AssignClassical('eCoordinateSystemJ2000',a,e,i,w,Raan,M);
propagator.Propagate; 

% ��������������Ӵ�����������������ÿ�������϶����д�����
seedSensor = seedSat.Children.New('eSensor','Sen'); 
% ���ô���������
seedSensor.CommonTasks.SetPatternSimpleConic(44.85, 0.1);

% ��������
disp('����walker����');
nPlan = 2;% ƽ����
nPerPlan = 2;% ÿ��ƽ��������
nRANNSpreed = 1;% ����ƽ��������λ��
% ����ѭ����Ÿ�ʽ���ã���ʾ��λ�����
nFormatPlan = 1;
nFormatPerPlan = 1;
strFormatPlan = ['%0' int2str(nFormatPlan) 'd'];
strFormatPerPlan = ['%0' int2str(nFormatPerPlan) 'd'];
% �Ҳ������󷽷��ˣ���Connect����
strWalkerSet = ['Walker */Satellite/Sat Delta ' int2str(nPlan) ' ' int2str(nPerPlan) ' ' int2str(nRANNSpreed) ' 360.0 No'];
root.ExecuteCommand(strWalkerSet);
% Ĭ�ϵ���������ֻ�����ǣ��������´���������������������д�����
constellation = root.CurrentScenario.Children.New('eConstellation','MyConst');
% ѭ����Ӵ�����
for i = 1:nPlan
    for j = 1:nPerPlan
        strSensorPath = ['*/Satellite/Sat' num2str(i,strFormatPlan) num2str(j,strFormatPerPlan)  '/Sensor/Sen'];
        constellation.Objects.Add(strSensorPath);
    end
end

% ������������
disp('������������');
covDef = scenario.Children.New('eCoverageDefinition','CovDef');
covDef.AssetList.Add(constellation.Path);
% ��������
disp('��������');
covDef.ComputeAccesses();

% ����Ʒ�ʲ���
disp('����FOM����');
fom = covDef.Children.New('eFigureofMerit','Fom');
fom.SetDefinitionType('eFmNAssetCoverage');
fom.Definition.Satisfaction.EnableSatisfaction = true;
% Find min/max FOM value for static contours
overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
Result_1 = overallValDP.Exec();
min = cell2mat(Result_1.DataSets.GetDataSetByName('Minimum').GetValues);
max = cell2mat(Result_1.DataSets.GetDataSetByName('Maximum').GetValues);

% end