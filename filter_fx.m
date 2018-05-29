function [ denosied_fx ] = filter_fx( org )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
lf =8;        % Operator length
flow = 1;       % Min freq to process
fhigh=125;      % Max freq to process180
mu = 0.9;      % trade-off parameter0.01
dt=0.004;

[nPoint,nShot]=size(org);

nBlockT =20; %ʱ����
nBlockX =20; %���᣻
dt=0.004;

nOverlapT = floor(nBlockT/2);
nOverlapX = floor(nBlockX/2);
% ������
% ˵�����Էֿ���ݣ����ú��������м�Ȩ��?
weightX = hanning(nBlockT);
weightY = hanning(nBlockX);
Weight2D = kron(weightX, weightY');
threshold_mul=1e-100;

% ����ֿ飬��ʱ����͵����ϵķֿ���
% ���У�nBx��ʱ�����ϵķֿ���nBy�ǵ����ϵķֿ���
nBx = ceil((nPoint - nOverlapT) / (nBlockT - nOverlapT));
nBy = ceil((nShot - nOverlapX) / (nBlockX - nOverlapX));

ePriL2 = zeros(nPoint,nShot);
ePriWeightL2 = zeros(nPoint,nShot);
    for iBx = 1 : nBx
        % ����ֿ���ʱ���᷽���ϵ���ʼλ��?
        xStart = (iBx - 1) * (nBlockT - nOverlapT) + 1;
        % ����ֿ���ʱ���᷽���ϵ���ֹλ��?
        xEnd = xStart + nBlockT - 1;
        % �����ֹλ�÷���Խ�磬����������ʱ��㣬
        % �����÷ֿ����ֹλ��������ʱ���
        % ����Ӧ�ĸı�ֿ����ʼλ��
        if xEnd > nPoint
            xEnd = nPoint;
            xStart = xEnd - nBlockT + 1;
        end

        for iBy = 1 : nBy
            % ����ֿ��ڵ��᷽���ϵ���ʼλ��?
            yStart = (iBy - 1) * (nBlockX - nOverlapX) + 1;
            % ����ֿ��ڵ��᷽���ϵ���ֹλ��?
            yEnd = yStart + nBlockX - 1;
            % �����ֹλ�÷���Խ�磬���������ĵ���?
            % �����÷ֿ����ֹλ����������?
            % ����Ӧ�ĸı�ֿ����ʼλ��
            if yEnd > nShot
                yEnd = nShot;
                yStart = yEnd - nBlockX + 1;
            end
            orgBlock = org(xStart:xEnd, yStart:yEnd);
            energyBlock = sum(sum(abs(orgBlock))) / (nBlockT * nBlockX);
            if energyBlock < threshold_mul
                ePriL2(xStart: xEnd, yStart : yEnd) = ePriL2(xStart: xEnd, yStart : yEnd) + orgBlock.* Weight2D;
                ePriWeightL2(xStart: xEnd, yStart : yEnd) = ePriWeightL2(xStart: xEnd, yStart : yEnd) + Weight2D;
                continue;
            end
    %         block_tv=tv_L2_wht(orgBlock,1,5);

            Filter_fx = fx_decon(orgBlock,dt,lf,mu,flow,fhigh);

            ePriL2(xStart: xEnd, yStart : yEnd) = ePriL2(xStart: xEnd, yStart : yEnd) + Filter_fx.* Weight2D;
            ePriWeightL2(xStart: xEnd, yStart : yEnd) = ePriWeightL2(xStart: xEnd, yStart : yEnd) + Weight2D;
        end 
    end
    denosied_fx = ePriL2 ./ ePriWeightL2;

end

