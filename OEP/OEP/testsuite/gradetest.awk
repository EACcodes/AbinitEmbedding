function abs(value)
{
  return (value<0?-value:value);
}

{ 
    if (abs($1 - $2) < 1e-5)
    {
	print $1,$2,"ok";
    } else { 
	print $1,$2,"fail"
    }
}