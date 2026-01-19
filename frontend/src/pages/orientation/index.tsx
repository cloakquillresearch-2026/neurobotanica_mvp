import { useEffect } from 'react';
import { useRouter } from 'next/router';

export default function OrientationIndex() {
  const router = useRouter();

  useEffect(() => {
    router.replace('/orientation/welcome');
  }, [router]);

  return null;
}
