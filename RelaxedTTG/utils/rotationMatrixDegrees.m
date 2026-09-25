function R = rotationMatrixDegrees(angle)
%ROTATIONMATRIXDEGREES Counterclockwise two-dimensional rotation.
c = cosd(angle);
s = sind(angle);
R = [c,-s;s,c];
end
