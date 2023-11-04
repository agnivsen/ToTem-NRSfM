N = 50;
nF = 5;

figure;
for ii = 1:nF
    DataCell{ii} = rand(3,N);
    DataCell{ii}(3,:) = zeros(1,N);
    VisVecCell{ii} = ones(1,N);
    plot3(DataCell{ii}(1,:),DataCell{ii}(2,:),DataCell{ii}(3,:),'k.'); hold on;
end

GPA = KernelGPA;%DefGPA('AFFINE');
GPA.centerDataAhead = true;
GPA.addDataByArray(DataCell, VisVecCell);
GPA.run();
GPA.mShape

trns = GPA.transformPoints(DataCell{3}, 3)

plot3(trns(1,:),trns(2,:),trns(3,:),'ro'); hold on;
plot3(GPA.mShape(1,:),GPA.mShape(2,:),GPA.mShape(3,:),'bo'); hold on;