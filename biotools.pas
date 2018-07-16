unit biotools;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, FileUtil, Forms, Controls, Graphics, Dialogs, Math,GraphType, ExtCtrls;
type
  Tpunto = record
    X,Y,Z:real;
  end;

  TPerfil= array of real;

  TTabla= array['A'..'Z'] of real;

  TAtomPDB = record
    NumAtom : integer;
    ID : string;
    AA : string;
    AltRes: char;
    Sub : char;
    residuo : integer;
    coor : Tpunto;
    Temp : Real;
  end;
  TInt=array of integer;
  Tperf=array of real;
  TPDB = record
    atm: array of TatomPDB;
    Header: string;
    totalatm: integer;
    totalres: integer;
    resindex: array of integer;
    end;
  Ttabladatos = array of array of real;
  Tpuente=array of integer;
  const
    borde=10;
    AA= 'ALA=A#CYS=C#ASP=D#GLU=E#PHE=F#GLY=G#HIS=H#ILE=I#LYS=K#LEU=L#MET=M#ASN=N#PRO=P#GLN=Q#ARG=R#SER=S#THR=T#VAL=V#TRP=W#TYR=Y';

function distance(p1,p2:Tpunto):real;
function loadPDB (txt:Tstrings) : TPDB;
function torsionAngle (p1,p2,p3,p4:Tpunto):real;
function bondAngle (p1,p2,p3:Tpunto):real;
function vectorialProduct (p1,p2:Tpunto):Tpunto;
function module (p:Tpunto):real;
function diffVector (p1,p2:Tpunto):Tpunto;
function sumVector(p1,p2:Tpunto):Tpunto;
function scalarProduct (p1,p2:Tpunto):real;
function ramachPlot (p:TPDB; var datos:Ttabladatos):boolean;
function atomSearch (var n1: integer; var n2: integer; p1:TPDB) : boolean;
function atomSearch (var n1: integer; var n2: integer; var n3:integer; p1:TPDB) : boolean; overload;
function atomSearch (var n1: integer; var n2: integer; var n3:integer; var n4:integer; p1:TPDB) : boolean;   overload;
function plot(datos: TTablaDatos; graf: TImage; cx: integer=1; cy:integer=2;
                clean:boolean=true; clpluma:TColor=clyellow;
                clrelleno:TColor=clyellow; clfondo:TColor=clblack;
                estilo: integer=0): boolean;
function rotate(p1:Tpunto;angulo:float;Eje:char):Tpunto;
function translate(p1:Tpunto;deltax:float;deltay:float;deltaz:float):Tpunto;
function rightscale(cadena:string;ancho:integer):string;
function leftscale(cadena:string;ancho:integer):string;
function savePDB(p:TPDB):Tstrings;
function mutation(p2:TPDB):TPDB;
function profile (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
function stroud (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
function eisenberg (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
function shortalong(cod3: char): string;
function longashort (cod3:string):char;
function cysteines(p1:TPDB): Tint;
function bridges(p1:TPDB;dist:real):Tstrings;
function RMSD(CysIndex:Tint;p1:TPDB;n1:integer;n2:integer):real;
function pheconectivity(p1:TPDB):TPDB;
implementation
function distance(p1,p2:Tpunto):real;
var
resultado:real;
begin
     resultado:=sqrt(sqr(p1.X-p2.X)+sqr(p1.Y-p2.Y)+sqr(p1.Z-p2.Z));
     distance:=resultado;
end;


function loadPDB (txt:Tstrings) : TPDB;
var
  p : TPDB;
  j,n,m,i,a,b : integer;
  linea : string;
  atomOK: boolean;

begin
  j:=0;
  n:=0;
  m:=0;
  a:=1;
  b:=1;

  for i:=1 to txt.Count-1 do
      begin
           linea:=txt[i];
           if pos ('ATOM', linea)=1 then
           a:=a+1;
           if trim(copy(linea,13,4))='CA' then
            b:=b+1;
      end;
          setlength(p.atm,a);
          setlength(p.resindex,b);
       begin
            for j:= 1 to txt.Count-1 do
            begin
                 linea:=txt[j];
            if pos ('ATOM', linea)=1 then
            begin
            atomOK:= (linea[17]=' ') OR
            ((linea[17])<>' ') AND  ((p.atm[n].residuo<>strtoint(trim(copy(linea,23,4)))) OR (linea[17]= p.atm[n].AltRes));
            begin
            n:=n+1;
            p.atm[n].NumAtom:= strtoint(trim(copy (linea,7,5)));
            p.atm[n].ID:= trim(copy(linea,13,4));
            p.atm[n].AltRes:= linea[17];
            p.atm[n].AA:=copy(linea,18,3);
            p.atm[n].Sub:= linea[22];
            p.atm[n].residuo:= strtoint(trim(copy(linea,23,4)));
            p.atm[n].coor.X := strtofloat(trim(copy(linea,31,8)));
            p.atm[n].coor.Y := strtofloat(trim(copy(linea,39,8)));
            p.atm[n].coor.Z := strtofloat(trim(copy(linea,47,8)));
            p.atm[n].temp   := strtofloat(trim(copy(linea,61,6)));

                              if p.atm[n].ID='CA' then
                              begin
                                   m:=m+1;
                                   p.resindex[m]:=n;
                              end;
            end
            end;
            end;
       end;

       p.totalatm:=n;
       p.totalres:=m;
       loadPDB:=p;
     end;
function module (p:Tpunto):real;
begin
    module:= sqrt(p.X*p.X + p.Y*p.Y + p.Z*p.Z);
end;
function diffVector (p1,p2:Tpunto):Tpunto;
var
  resultado: Tpunto;
begin
    resultado.X:=p2.X-p1.X;
    resultado.Y:=p2.Y-p1.Y;
    resultado.Z:=p2.Z-p1.Z;
    diffVector:=resultado;
end;
function sumVector(p1,p2:Tpunto):Tpunto;
var
  suma:Tpunto;
begin
  suma.x:=p2.x+p1.x;
  suma.y:=p2.y+p1.y;
  suma.z:=p2.z+p1.z;
  sumVector:=suma;
end;
function scalarProduct (p1,p2:Tpunto):real;
begin
    scalarProduct:=p1.X*p2.X+p1.Y*p2.Y+p1.Z*p2.Z;
end;
   function bondAngle (p1,p2,p3:Tpunto):real;
   var
     v1,v2:Tpunto;
   begin
       v1:=vectorDif(p2,p1);
       v2:=vectorDif(p2,p3);
       bondAngle:=radtodeg(arccos(prodEscalar(v1,v2)/modulo(v1)*modulo(v2)));

   end;
   function vectorialProduct (p1,p2:Tpunto):Tpunto;
   var
     resultado: Tpunto;
     begin
         resultado.X:=p1.Y*p2.Z-p1.Z*p2.Y;
         resultado.Y:=p1.Z*p2.X-p1.X*p2.Z;
         resultado.Z:=p1.X*p2.Y-p1.Y*p2.X;
         vectorialProduct:=resultado;
     end;
   function torsionAngle (p1,p2,p3,p4:Tpunto):real;
   var
     v1,v2,v,BA,BC,CB,CD:Tpunto;
     ang,prod,arco: real;

     begin
         BA:= vectorDif(p1,p2);
         BC:= vectorDif(p2,p3);
         {CB:= vectorDif(p3,p2); }
         CD:= vectorDif(p3,p4);
         v1:=prodVectorial(BA,BC);
         v2:=prodVectorial(BC,CD);
         v:=prodVectorial(v2,v1);

         arco:= arccos(prodEscalar(v1,v2)/(modulo(v1)*modulo(v2)));
         prod:= prodEscalar(BC,v)/(modulo(BC)*modulo(v));
         if prod>0 then prod:=1 else prod:= -1;
         if abs(prodEscalar(v1,v2)/(modulo(v1)*modulo(v2)))>1 then showmessage('Coseno fuera del rango real') else
         ang:=arco*prod;
         torsionAngle:= radtodeg(ang);
     end;
   function atomSearch (var n1: integer; var n2: integer; p1:TPDB) : boolean;
   var
      j: integer;
     flag1,flag2: boolean;
   begin
     flag1:=true; flag2:=true;
     if (p1.atm[n1].NumAtom<>n1) OR (p1.atm[n2].NumAtom<>n2)  then
     begin
       flag1:=false; flag2:=false;
       for j:=1 to p1.totalatm do
       begin
         if p1.atm[j].NumAtom=n2 then
         begin
              n2:=j;
              flag2:=true;
         end;
       end;
       begin
         if p1.atm[j].NumAtom=n1 then
         begin
              n1:=j;
              flag1:=true;
         end;
       end;
       end;
   if flag1 AND flag2 then buscarATMnum:= true
      else
        begin
          showmessage ('Alguno de los átomos no se encuentra. Vuelva a introducirlos');
          atomSearch:= false;
        end;
   end;
   function atomSearch (var n1: integer; var n2: integer; var n3: integer; p1:TPDB) : boolean;
   var
      j: integer;
     flag1,flag2,flag3: boolean;
   begin
        flag1:=true; flag2:=true; flag3:=true;
        if (p1.atm[n1].NumAtom<>n1) OR (p1.atm[n2].NumAtom<>n2)OR (p1.atm[n3].NumAtom<>n3)  then
        begin
        flag1:=false; flag2:=false; flag3:=false;
        for j:=1 to p1.totalatm do
        begin
             if p1.atm[j].NumAtom=n2 then
             begin
                  n2:=j;
                  flag2:=true;
             end;
        end;
        begin
             if p1.atm[j].NumAtom=n1 then
             begin
                  n1:=j;
                  flag1:=true;
             end;
        end;
        begin
             if p1.atm[j].NumAtom=n3 then
             begin
                  n3:=j;
                  flag3:=true;
             end;
        end;
        end;
   if flag1 AND flag2 AND flag3 then buscarATMnum:= true
      else
        begin
          showmessage ('Alguno de los átomos no se encuentra. Vuelva a introducirlos');
          atomSearch:= false;
        end;
   end;
      function atomSearch (var n1: integer; var n2: integer; var n3: integer; var n4: integer; p1:TPDB) : boolean;
   var
      j: integer;
     flag1,flag2,flag3, flag4: boolean;

   begin
        flag1:=true; flag2:=true; flag3:=true; flag4:=true;
        if (p1.atm[n1].NumAtom<>n1) OR (p1.atm[n2].NumAtom<>n2)OR (p1.atm[n3].NumAtom<>n3) OR (p1.atm[n4].NumAtom<>n4) then
        begin
            flag1:=false; flag2:=false; flag3:=false; flag4:=false;
            for j:=1 to p1.totalatm do
            begin
                 if p1.atm[j].NumAtom=n2 then
                 begin
                      n2:=j;
                      flag2:=true;
                 end;
            end;
            begin
                 if p1.atm[j].NumAtom=n1 then
                 begin
                      n1:=j;
                      flag1:=true;
                 end;
            end;
            begin
                 if p1.atm[j].NumAtom=n3 then
                 begin
                      n3:=j;
                      flag3:=true;
                 end;
            end;
            begin
                 if p1.atm[j].NumAtom=n4 then
                 begin
                      n4:=j;
                      flag4:=true;
                 end;
            end;
        end;
   if flag1 AND flag2 AND flag3 AND flag4 then buscarATMnum:= true
   else
        begin
          showmessage ('Missing atoms');
          atomSearch:= false;
        end;
   end;
function ramachPlot (p:TPDB; var datos:Ttabladatos):boolean;
var
  CP,N,CA,C,NP: Tpunto;
  j:integer;
  OK: boolean;

  begin
      OK:= true;
      setlength(datos,3,p.totalres);
        for j:=2 to p.totalres-1 do
          begin
            CA:=p.atm[p.resindex[j]].coor; if p.atm[p.resindex[j]].ID<> 'CA' then OK:= false;
            N:=p.atm[p.resindex[j]-1].coor; if p.atm[p.resindex[j]-1].ID <> 'N' then OK:= false;
            C:=p.atm[p.resindex[j]+1].coor; if p.atm[p.resindex[j]+1].ID <> 'C' then OK:= false;
            CP:=p.atm[p.resindex[j-1]+1].coor; if p.atm[p.resindex[j-1]+1].ID <> 'C' then ok:=false;
            NP:=p.atm[p.resindex[j+1]-1].coor; if p.atm[p.resindex[j+1]-1].ID <>'N' then OK:=false;
            if OK then
            begin
              datos [1,j-2]:= angulotorsion (CP,N,CA,C);
              datos [2,j-2]:= angulotorsion (N,CA,C,NP);
            end;
          end;
        ramachPlot := true;
   end;

function plot(datos: TTablaDatos; graf: TImage; cx: integer=1; cy:integer=2;
                clean:boolean=true; clpluma:TColor=clyellow;
                clrelleno:TColor=clyellow; clfondo:TColor=clblack;
                estilo: integer=0): boolean;

var
  xmin, ymin: real;
  rangox, rangoy: real;
  j, xx, yy, ancho, alto: integer;
  tipolinea:TPenStyle;

  function fitx(x:real): real;
  begin
    fitx:= (ancho-2*borde)*(x-xmin)/(rangox);
  end;

  function fity(y:real): real;
  begin
    fity:= (alto-2*borde) - ((alto-2*borde)*(y-ymin)/(rangoy));
  end;

begin
  if (high(datos)>=cx) and (high(datos)>=cy) and (high(datos[1])>=0) then
  begin
    xmin:= minvalue(datos[cx]);
    ymin:= minvalue(datos[cy]);
    rangox:= maxvalue(datos[cx])-xmin;
    rangoy:= maxvalue(datos[cy])-ymin;
    ancho:= graf.Width-1; alto:= graf.Height-1;
    graf.Canvas.Pen.Style:=psSolid;
    if clean then
    begin
      graf.canvas.Pen.Color:= clfondo;
      graf.canvas.Brush.Color:= clfondo;
      graf.canvas.Rectangle(0,0,graf.Width-1, graf.Height-1);
    end;
    graf.canvas.Pen.Color:= clpluma;
    graf.canvas.Brush.Color:= clrelleno;

   for j:=0 to high(datos[0]) do
      begin
        xx:= round(fitx(datos[cx,j]));
        yy:= round(fity(datos[cy,j]));
        graf.canvas.Ellipse(xx-2, yy-2, xx+2, yy+2);
      end;
    if estilo=1 then tipolinea:= psSolid else tipolinea:= psClear;
    graf.Canvas.Pen.Style:=tipolinea;
    graf.Canvas.MoveTo(round(fitx(datos[cx,0])), round(fity(datos[cy,0])));
    for j:=0 to high(datos[0]) do
      begin
        xx:= round(fitx(datos[cx,j]));
        yy:= round(fity(datos[cy,j]));
        graf.canvas.LineTo(xx, yy);
      end;
    plot:= true;
    end
    else
    begin
      plot:= false;
      showmessage('Error plot');
    end;
end;
function translate(p1:Tpunto;deltax:float;deltay:float;deltaz:float):Tpunto;
var
  puntonuevo:Tpunto;
begin
  puntonuevo.x:=p1.x+deltax;
  puntonuevo.y:=p1.y+deltay;
  puntonuevo.z:=p1.z+deltaz;
  translate:=puntonuevo;
end;

function Rotacion(p1:Tpunto;angulo:float;Eje:char):Tpunto;
var
  puntonuevo:Tpunto;
begin
 if Eje='x' then
 begin
   puntonuevo.x:=p1.x;
   puntonuevo.y:=p1.y*cos(angulo)-p1.z*sin(angulo);
   puntonuevo.z:=p1.y*sin(angulo)+p1.z*cos(angulo);
 end;
 if Eje='y' then
 begin
   puntonuevo.x:=p1.x*cos(angulo)+p1.z*sin(angulo);
   puntonuevo.y:=p1.y;
   puntonuevo.z:=p1.z*cos(angulo)-p1.x*sin(angulo);
 end;
 if Eje='z' then
 begin
   puntonuevo.x:=p1.x*cos(angulo)-p1.y*sin(angulo);
   puntonuevo.y:=p1.x*sin(angulo)+p1.y*cos(angulo);
   puntonuevo.z:=p1.z;
 end;
 rotation:=puntonuevo;
end;

function righscale(cadena:string;ancho:integer):string;
var
  resultado:string;
begin
 resultado:='                                                   '+cadena;
 resultado:=copy(resultado,length(resultado)-ancho+1,ancho);
 rightscale:=resultado;
end;

function leftscale(cadena:string;ancho:integer):string;
var
  resultado:string;
begin
 resultado:=cadena+'                                                   ';
 resultado:=copy(resultado,1,ancho);
 leftscale:=resultado;
end;

function savePDB(p:TPDB):Tstrings;
var
  txt:Tstrings;
  j:integer;
  linea:string;
begin
 txt:=TstringList.Create;
 for j:=1 to p.totalatm-1 do
 begin
 linea:='ATOM  '+escaladerecha(inttostr(p.atm[j].NumAtom),5)+'  '+escalaizquierda(p.atm[j].ID,3)+p.atm[j].AltRes+p.atm[j].AA+
   ' '+p.atm[j].sub+escaladerecha(inttostr(p.atm[j].residuo),4)+'    '+escaladerecha(formatfloat('0.000',p.atm[j].coor.x),8)+
   escaladerecha(formatfloat('0.000',p.atm[j].coor.y),8)+escaladerecha(formatfloat('0.000',p.atm[j].coor.z),8)+ escaladerecha(formatfloat('0.000',p.atm[j].Temp),8);
   txt.Add(linea);
 end;
savePDB:=txt;
end;
function shortalong(cod3: char): string;
var
  posicion: integer;
  residuo: string;
  begin
  posicion:= pos(uppercase(cod3),AA);
  if posicion > 0 then residuo:= copy(AA,(posicion+2),3) else residuo:='***';
  shortalong:= residuo
  end;

function longashort (cod3:string):char;
var
  {posicion: integer;}
  letra: char;
  begin
  cod3:=uppercase(cod3);
  if cod3='HIS' then letra:='H'

  else if cod3='ARG' then letra:='R'
  else if cod3='LYS' then letra:='K'
  else if cod3='ILE' then letra:='I'
  else if cod3='PHE' then letra:='F'
  else if cod3='LEU' then letra:='L'
  else if cod3='TRP' then letra:='W'
  else if cod3='ALA' then letra:='A'
  else if cod3='MET' then letra:='M'
  else if cod3='PRO' then letra:='P'
  else if cod3='CYS' then letra:='C'
  else if cod3='ASN' then letra:='N'
  else if cod3='VAL' then letra:='V'
  else if cod3='GLY' then letra:='G'
  else if cod3='SER' then letra:='S'
  else if cod3='GLN' then letra:='Q'
  else if cod3='TYR' then letra:='Y'
  else if cod3='ASP' then letra:='D'
  else if cod3='GLU' then letra:='E'
  else if cod3='THR' then letra:='T'
  else if cod3='0' then cod3:=''
  else if cod3='?' then cod3:=''
  else
  letra:='*';
  {posicion:= pos(uppercase(cod3),AA);
  if posicion > 0 then letra  := AA[posicion-2] else letra:='*'; }
  longashort:= letra;
  end;
function cysteines(p1:TPDB):TInt;
var
 Indice:TInt;
 j,i,n,m:integer;
begin

 n:=0;
 for j:=1 to p1.totalatm do
 begin
  if p1.atm[j].ID='SG' then n:=n+1;
 end;
 setlength(Indice,n+1);

 m:=1;
 for i:=1 to p1.totalatm do
 begin
  if p1.atm[i].ID='SG' then
  begin
   Indice[m]:=i;
   m:=m+1;
  end;
 end;

cysteines:=Indice;
end;
function RMSD(CysIndex:TInt;p1:TPDB;n1:integer;n2:integer):real;
var
 Cys:Tperfil;
 resultado:real;
 j:integer;
begin
 setlength(Cys,7);
 for j:=1 to 6 do Cys[j]:=sqr(dist3D(p1.atm[CysIndex[n1]-(j-1)].coor,p1.atm[CysIndex[n2]-(j-1)].coor));
 resultado:=sqrt((Cys[1]+Cys[2]+Cys[3]+Cys[4]+Cys[5]+Cys[6])/6);
 RMSD:=resultado;
end;
function mutation(p2:TPDB):TPDB;
var
 p1:TPDB;
 n,j,i,k,total:integer;
 v1,v2,v3,v4,v5,vu,vf:Tpunto;
 m3:real;
begin

 p1:=p2;
 total:=p1.totalatm;
 n:=1;
 setlength(p1.atm,p1.totalatm+2);
 while p1.atm[n].AA <> 'PHE' do n:=n+1;
 while p1.atm[n].AA = 'PHE' do n:=n+1;

 for j:=total downto n do p1.atm[j+1]:=p1.atm[j];

 for i:=n downto n-11 do p1.atm[i].AA:='TYR';

 p1.atm[n].ID:='OH';

 for k:=total+1 downto n+1 do p1.atm[k].NumAtom:=p1.atm[k].NumAtom+1;

 p1.atm[n].residuo:=p1.atm[n-1].residuo;

 v1:=p1.atm[n-7].coor;
 v2:=p1.atm[n-6].coor;
 v5:=p1.atm[n-1].coor;
 v3:=VectorDif(v1,v2);
 m3:=Modulo(v3);
 vu.x:=v3.x/m3;
 vu.y:=v3.y/m3;
 vu.z:=v3.z/m3;
 v4.x:=vu.x*1.38;
 v4.y:=vu.y*1.38;
 v4.z:=vu.z*1.38;
 vf:=VecSum(v4,v5);
 p1.atm[n].coor:=vf;

 mutation:=p1;

end;
function profile (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
var
   j,k, ndatos : integer;
   resultado: TPerfil;
   suma, media: real;
begin
     setlength(resultado,length(seq)+1);
      for j:=1 to length(seq) do
      begin
           suma:= 0;
           ndatos:= 0;
           for k:= max(j-semiventana,1) to min(j+semiventana,length(seq)) do
           begin
                suma:= suma+tabla[seq[k]];
                ndatos:= ndatos+1;
           end;
           media:= suma/ndatos;
           resultado[j]:=media;
      end;
profile:=resultado;
end;
function stroud (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
var
   j,k,i, ndatos : integer;
   resultado: TPerfil;
   suma,suma1, suma2, delta, media, final: real;
begin
     setlength(resultado,length(seq)+1);
      for j:=1 to length(seq) do
      begin
           delta:=degtorad(100);
           suma:= 0;
           suma1:= 0;
           suma2:= 0;
           ndatos:= 0;
           for k:= max(j-semiventana,1) to min(j+semiventana,length(seq)) do
           begin
                suma:= suma+tabla[seq[k]];
                ndatos:= ndatos+1;
           end;
           media:= suma/ndatos;
           for i:=max(j-semiventana,1) to min(j+semiventana,length(seq)) do
           begin
           suma1:=suma1+((tabla[seq[i]]-media)*sin(delta*i));
           suma2:=suma2+((tabla[seq[i]]-media)*cos(delta*i));
           end;
           final:=sqr(suma1)+sqr(suma2);
           resultado[j]:=final;
      end;
     stroud:=resultado;
end;
function eisenberg (seq: string; tabla: TTabla; semiventana: integer): TPerfil;
var
   j,k, ndatos : integer;
   resultado: TPerfil;
   suma1, suma2, media, delta: real;
begin
     setlength(resultado,length(seq)+1);
      for j:=1 to length(seq) do
      begin
           suma1:= 0;
           suma2:= 0;
           ndatos:= 0;
           delta:= degtorad(100);
           for k:= max(j-semiventana,1) to min(j+semiventana,length(seq)) do
           begin
                suma1:= suma1+tabla[seq[k]]*sin(delta*k);
                suma2:= suma2+tabla[seq[k]]*cos(delta*k);
                ndatos:= ndatos+1;
           end;
             media:=sqrt(sqr(suma1)+sqr(suma2));
             resultado[j]:=media;
      end;
eisenberg:=resultado;
end;
function pheconectivity(p1:TPDB):TPDB;
var
 prot:TPDB;
 j,n,i,m:integer;
begin
setlength(prot.atm,500);
j:=1;
n:=1;

 while p1.atm[j].AA <> 'PHE' do j:=j+1;
 while p1.atm[j].AA = 'PHE' do
 begin
  if p1.atm[j].ID = 'CG' then m:=j;
  prot.atm[n]:=p1.atm[j];
  n:=n+1;
  j:=j+1;
 end;

 for i:=1 to p1.totalatm-1 do
 begin
  if (Dist3D(p1.atm[m].coor,p1.atm[i].coor)<10) and (Dist3D(p1.atm[m].coor,p1.atm[i].coor)>0) then
  begin
   prot.atm[n]:=p1.atm[i];
   n:=n+1;
  end;
 end;
 prot.totalatm:=n;
 pheconectivity:=prot;
end;
function bridges(p1:TPDB;dist:real):Tstrings;
var
 Lista: TInt;
 Matriz: Ttabladatos;
 i,j, k, n : integer;
 Linea:string;
 txt:Tstrings;
begin
 Lista:=TablaCys(p1);
 setlength(Matriz,length(Lista)+1,length(Lista)+1);
 for j:=1 to length(Lista)-1 do
 begin
  for i:=1 to length(Lista)-1 do
     Matriz[j,i]:=Dist3D(p1.atm[Lista[j]].coor,p1.atm[Lista[i]].coor);
 end;
 begin
  txt:=Tstringlist.create;
  for k:=1 to length(Lista) do
  begin
       for n:=k+1 to length(Lista)-1 do
       begin
            if Matriz[k,n]<dist then
            begin
                 Linea:='Disulfide bridge between: cysteine '+inttostr(p1.atm[Lista[k]].residuo)+' and '+'cysteine '+inttostr(p1.atm[Lista[n]].residuo);
                 txt.add(Linea);
            end;
            bridges:= txt;
       end;
  end;
 end;
end;
end.
