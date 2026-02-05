import type { ConditionPayload } from '@/types/customer'

export type ConditionInput = string | ConditionPayload | null | undefined

export function conditionName(condition: ConditionInput): string {
  if (!condition) return ''
  if (typeof condition === 'string') return condition
  return condition.name || ''
}

export function conditionNames(conditions?: ConditionInput[]): string[] {
  return (conditions ?? [])
    .map(conditionName)
    .filter((name): name is string => Boolean(name && name.length))
}
