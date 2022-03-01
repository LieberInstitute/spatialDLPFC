function [O] = extractMD(fname)

imgInfo = imfinfo(fname);
N = size(imgInfo,1);

for i = 1:N-1
temp = imgInfo(i).ImageDescription;
% temp = strrep(temp,'utf-16','utf-8');
% xmlDocument = javax.xml.parsers.DocumentBuilderFactory.newInstance().newDocumentBuilder.parse(java.io.StringBufferInputStream(temp));
% nodeList = xmlDocument.getElementsByTagName('Name');
% numberOfNodes = nodeList.getLength();
% firstNode = nodeList.item(0);
% O{i} = char(firstNode.getTextContent);
% O{i} = strtok(O{i},' ');
% end

T=regexp(temp,'<Name>', 'split');
T2=regexp(T{2},'</Name>', 'split');
O{i} = strtok(T2{1},' ');
end