% 调用stk分析Starlink星座
% 戴正旭 2019.06.23
% function StarlinkAnalyze_COM
% 打开STK，不新建场景
% Requires STK, STK Coverage and STK Integration licenses.
clc
clear all
% % 
% % % Initialize 
disp('打开STK程序，创建一个场景');
try
    % 获取运行中的STK实例句柄
    uiapp = actxGetRunningServer('STK11.application');
    root = uiapp.Personality2;
    checkempty = root.Children.Count;
    if checkempty == 0
        % 如果未发现场景，新建一个
        uiapp.visible = 1;
        root.NewScenario('StarLink_COM');
        scenario = root.CurrentScenario;
    else
        % 如果发现打开着的场景，询问是否关闭重建
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
    % STK没有启动，新建实例获取句柄
    uiapp = actxserver('STK11.application');
    root = uiapp.Personality2;
    uiapp.visible = 1;
    root.NewScenario('StarLink_COM');
    scenario = root.CurrentScenario;
end

% 设定场景时间
disp('设定场景时间');
strBegTime = '27 May 2019 06:14:00.000';
strEndTime = '28 May 2019 06:14:00.000';
scenario.SetTimePeriod(strBegTime,strEndTime);
% 设定场景动画开始时间
scenario.Epoch = strBegTime;
% 设定动画时间回到起始点，这里用的是Connect命令
root.ExecuteCommand('Animate * Reset');

% 新建种子卫星
disp('创建种子卫星');
seedSat = scenario.Children.New('eSatellite', 'Sat');
seedSat.SetPropagatorType('ePropagatorJ4Perturbation');

% 设置卫星根数
% 注意这里单位与mexConnect中不同，半长轴是公里,角度是度
a=6928.137;  e=0.0; i=53.0;w=0; Raan=160; M=0;
propagator = seedSat.Propagator; 
propagator.StartTime = scenario.StartTime;  % Set to scenario start time
propagator.StopTime = scenario.StopTime;    % Set to scenario stop time
propagator.Step = 60.0;
propagator.InitialState.Representation.AssignClassical('eCoordinateSystemJ2000',a,e,i,w,Raan,M);
propagator.Propagate; 

% 在种子卫星上添加传感器，创建星座后每个卫星上都会有传感器
seedSensor = seedSat.Children.New('eSensor','Sen'); 
% 设置传感器参数
seedSensor.CommonTasks.SetPatternSimpleConic(44.85, 0.1);

% 创建星座
disp('生成walker星座');
nPlan = 2;% 平面数
nPerPlan = 2;% 每个平面卫星数
nRANNSpreed = 1;% 相邻平面卫星相位差
% 后续循环编号格式化用，表示几位数填充
nFormatPlan = 1;
nFormatPerPlan = 1;
strFormatPlan = ['%0' int2str(nFormatPlan) 'd'];
strFormatPerPlan = ['%0' int2str(nFormatPerPlan) 'd'];
% 找不到对象方法了，用Connect命令
strWalkerSet = ['Walker */Satellite/Sat Delta ' int2str(nPlan) ' ' int2str(nPerPlan) ' ' int2str(nRANNSpreed) ' 360.0 No'];
root.ExecuteCommand(strWalkerSet);
% 默认的星座还是只有卫星，这里重新创建星座集，用来存放所有传感器
constellation = root.CurrentScenario.Children.New('eConstellation','MyConst');
% 循环添加传感器
for i = 1:nPlan
    for j = 1:nPerPlan
        strSensorPath = ['*/Satellite/Sat' num2str(i,strFormatPlan) num2str(j,strFormatPerPlan)  '/Sensor/Sen'];
        constellation.Objects.Add(strSensorPath);
    end
end

% 创建覆盖区域
disp('创建覆盖区域');
covDef = scenario.Children.New('eCoverageDefinition','CovDef');
covDef.AssetList.Add(constellation.Path);
% 分析覆盖
disp('分析覆盖');
covDef.ComputeAccesses();

% 创建品质参数
disp('创建FOM参数');
fom = covDef.Children.New('eFigureofMerit','Fom');
fom.SetDefinitionType('eFmNAssetCoverage');
fom.Definition.Satisfaction.EnableSatisfaction = true;
% Find min/max FOM value for static contours
overallValDP = fom.DataProviders.GetDataPrvFixedFromPath('Overall Value');
Result_1 = overallValDP.Exec();
min = cell2mat(Result_1.DataSets.GetDataSetByName('Minimum').GetValues);
max = cell2mat(Result_1.DataSets.GetDataSetByName('Maximum').GetValues);

% end