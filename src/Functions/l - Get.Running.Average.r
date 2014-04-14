Get.Running.Average = function( v, n=length(v)/100 )
{
                               
	ret = c( NA, length(v) )


	window = sum ( v[ 1 : ( 2 * n + 1 ) ] )
	
	for( i in 1:length(v) )
	{
		if( i > n + 1 &&
		    i <= length(v) - n )
                           	{                                        
			window = window + v[ i + n ]
			window = window - v[ i - n - 1 ]
		}

		ret[ i ] = window / ( 2 * n + 1 )			

	}

	names( ret ) = names( v )
	return( ret )

}

